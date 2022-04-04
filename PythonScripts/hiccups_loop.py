import subprocess
import pandas as pd
import os

pairs = ['OM', 'OF', 'YM', 'YF']
hiccups_output = '../Results/loops/HICCUPS/{}'

for i in pairs:
    if not os.path.exists(hiccups_output.format(i)):
        print("Creating hiccups loops")
        subprocess.call('java -jar /data/pycharm/BML_HiC_Data_Analysis/juicer_tools_1.22.01.jar'
                        ' hiccups -m 500 -r 10000 '
                        '-f 0.1 -p 2 -i 5 -d 20000  '
                        '/data/MuSC_data/HiC_CPB_normalized/{}_new_CPBnorm_inter_30.hic  {}'.format(i, hiccups_output.format(i)),
                        shell=True)


exit()

sheet = 1
sheets = ['O_H3K122ac_enhancer_unique', 'O_H3K27ac_enhancer_unique']
enhancer_df = pd.read_excel("Unique_enhancer_list.xlsx", sheet_name=sheets[sheet])

regions = []
for line_number, (index, row) in enumerate(enhancer_df.iterrows()):
    start = row['START']
    stop = row['STOP']
    chrom = row['CHROM']
    regions.append((chrom, start, stop))
print(regions)

columns = ['Chrom', 'Enhancer_start', 'Enhancer_end', 'Loop_bin1_start',
           'Loop_bin1_end', 'Loop_bin2_start', 'Loop_bin2_end', 'type', 'FDR']
loops_file = ['output.loop1', 'output.loop2', 'output.diffloop1', 'output.diffloop2']
writer = pd.ExcelWriter('{}.xlsx'.format(sheets[sheet]), engine='xlsxwriter')

for loop_file in loops_file:
    print("Loops file is {}".format(loop_file))
    loop_df = pd.read_csv('/data/pycharm/BML_HiC_Data_Analysis/Results/loops/mustache_loops/{}'.format(loop_file)
                          , sep='\t')

    data = []
    for index, row in loop_df.iterrows():
        START_1 = row['BIN1_START']
        STOP_1 = row['BIN1_END']
        CHROM_1 = row['BIN1_CHR']
        START_2 = row['BIN2_START']
        STOP_2 = row['BIN2_END']
        CHROM_2 = row['BIN2_CHROMOSOME']
        FDR = row['FDR']

        assert (CHROM_1 == CHROM_2)

        # may be redundant, will look through
        for chrom, start, stop in regions:
            if chrom == CHROM_1:
                if start == START_1 and stop == STOP_1:
                    data.append([chrom, start, stop, START_1, STOP_1, START_2, STOP_2, 'same region', FDR])
                elif start == START_2 and stop == STOP_2:
                    data.append([chrom, start, stop, START_2, STOP_2, START_1, STOP_1, 'same region', FDR])

                elif start > START_1 and stop < STOP_1:
                    data.append([chrom, start, stop, START_1, STOP_1, START_2, STOP_2, 'inside region', FDR])
                elif start > START_2 and stop < STOP_2:
                    data.append([chrom, start, stop, START_2, STOP_2, START_1, STOP_1, 'inside region', FDR])

                elif start < START_1 and stop > STOP_1:
                    data.append([chrom, start, stop, START_1, STOP_1, START_2, STOP_2, 'outer region', FDR])
                elif start < START_2 and stop > STOP_2:
                    data.append([chrom, start, stop, START_2, STOP_2, START_1, STOP_1, 'outer region', FDR])

                elif start < START_1 < stop < STOP_1:
                    data.append([chrom, start, stop, START_1, STOP_1, START_2, STOP_2, 'left region', FDR])
                elif start < START_2 < stop < STOP_2:
                    data.append([chrom, start, stop, START_2, STOP_2, START_1, STOP_1, 'left region', FDR])

                elif START_1 < start < STOP_1 < stop:
                    data.append([chrom, start, stop, START_1, STOP_1, START_2, STOP_2, 'right region', FDR])
                elif START_2 < start < STOP_2 < stop:
                    data.append([chrom, start, stop, START_2, STOP_2, START_1, STOP_1, 'right region', FDR])

    _df = pd.DataFrame(columns=columns, data=data)
    _df.to_excel(writer, sheet_name=loop_file)

writer.save()
