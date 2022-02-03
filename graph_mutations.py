#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def graph(errors, total, output):
    xv = []
    yv = []
    for k,v in errors.items():
        xv.append(str(k) + '\n' + str(v[0]))
        yv.append(float(v[1]))
    plt.figure(figsize = [11,5])
    panel1 = plt.axes([1/10,.5/5,8/10,4/5])
    panel1.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
    sns.barplot(x=xv, y=yv, color = 'teal', ax = panel1)
    panel2 = plt.axes([7.8/10,3.8/5,1/10,.5/5])
    panel2.tick_params(axis = "both", which = "both", bottom = False, labelbottom = False, left = False, labelleft = False, right = False, labelright = False, top = False, labeltop = False)
    panel2.text(x = 1/10, y = 1.1/5, s ='Total Bases\n {}'.format(int(total)))
    plt.savefig(output, dpi = 600)

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help = 'path to input mutation file for graphing')
    parser.add_argument('-o', '--output', help = 'output name. default is input name with .png extension instead of .txt', default = None)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    errors = {}
    with open(args.input, 'r') as ein:
        for entry in ein:
            if entry[0] != "(":
                continue
            else:
                spent = entry.strip().split()
                key = spent[0][2] + '->' + spent[1][1]
                count = int(spent[2])
                rate = float(spent[3])
                errors[key] = (count, rate)
        #calculate genome size as well, using the last count/rate.
        total = sum([count/rate for count, rate in errors.values()])
    if args.output == None:
        output = args.input.split('.')[0] + '.png'
    elif args.output.endswith('.png'):
        output = args.output
    else:
        output = args.output + '.png'
    graph(errors, total, output)


if __name__ == "__main__":
    main()