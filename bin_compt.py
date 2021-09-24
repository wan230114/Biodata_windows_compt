#!/usr/bin/env python3
# -*- coding:utf-8 -*-

#############################################
# @ Author: Chen Jun
# @ Author Email: 1170101471@qq.com
# @ Created Date: 2021-09-22, 12:33:09
# @ Modified By: Chen Jun
# @ Last Modified: 2021-09-24, 22:34:06
#############################################

"""
1M为窗口，100kb为步长，进行滑动统计数据的峰值与平均值。
对每一个点计算上下窗口，统计数据情况，最后计算完一条染色体按染色体导出。
"""

# %%

infile = "./zma.index.chr1.chr2.xls"
outfile = "./zma.index.chr1.chr2.tongji.xls"

window_size = 10 ** 6  # 1M
step_size = 10 ** 5  # 100kb

fo = open(outfile, "w")


# %%
# 1-100  101-200  201-300
def get_windows(POS):
    L = []
    l_max = POS // step_size * step_size
    r_min = (POS // step_size + 1) * step_size
    # print("l_max:", l_max, "r_min:", r_min)
    LEFT, RIGHT = r_min-window_size+1, r_min
    if LEFT < 0:
        RIGHT = RIGHT + 1 - LEFT
        LEFT = 1
    while LEFT <= l_max:
        L.append((LEFT, RIGHT))
        LEFT += step_size
        RIGHT += step_size
    L.append((LEFT, RIGHT))
    return L


# test
get_windows(2000607)
# %%


def deal_D(D_chrs, old_CHR):
    tmp = D_chrs.pop(old_CHR)  # 删除元素并导出到文件
    for x in tmp:
        print(old_CHR, *x, *[resx for res in tmp[x]
                             for resx in res], sep="\t", file=fo)


old_CHR = ""
D_chrs = {old_CHR: {}}

with open(infile) as fi:
    fi.readline()
    i = 0
    for line in fi:
        Lline = line.strip().split("\t")
        # print(Lline)
        CHR, POS = Lline[0], int(Lline[1])
        print(CHR, POS)
        DATAs = [float(x) for x in Lline[2:]]
        if CHR not in D_chrs:
            print(f"dealing {CHR}")
            D_chr = {}
            D_chrs[CHR] = D_chr
            deal_D(D_chrs, old_CHR)
            old_CHR = CHR
        for loc in get_windows(POS):
            print("loc:", loc)
            if loc not in D_chr:
                MAX, MIN, SUM, NUM, AVG = 0, 0, 0, 0, 0
                D_chr[loc] = [[MAX, MIN, SUM, NUM, AVG]
                              for i in range(len(DATAs))]
            for x, comp in zip(DATAs, D_chr[loc]):
                print(x, comp)
                if x > comp[0]:
                    comp[0] = x
                if x < comp[0]:
                    comp[0] = x
                comp[2] += x
                comp[3] += 1
                comp[4] = comp[2]/comp[3]

deal_D(D_chrs, old_CHR)
