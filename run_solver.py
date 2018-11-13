import time
from matplotlib import pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import subprocess

from colorama import init, Fore, Back, Style
init()

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("command")
args = parser.parse_args()

cmd = args.command

t_init = 1.2

def find_between(s, first, last):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""

# draw the figure
step_list = [0]
t_max_list = [0]
q_he4_cu_list = [0]
q_cu_he3_list = [0]
plot_t_max, = plt.plot(step_list, t_max_list, 'red')
plot_q_he4_cu, = plt.plot(step_list, q_he4_cu_list, 'blue')
plot_q_cu_he3, = plt.plot(step_list, q_cu_he3_list, 'green')

plt.draw()
plt.ylim([0.00001, 10])
plt.yscale('log')
plt.pause(0.00001)

popen = subprocess.Popen(cmd, stdout=subprocess.PIPE)
lines_iterator = iter(popen.stdout.readline, b"")
i = 0
q_deposit = 0
while popen.poll() is None:
    for line in lines_iterator:
        nline = line.rstrip().decode("latin")
        if not q_deposit:
            q_deposit = find_between(nline, 'q_deposit = ', ';').strip(' ')
        step = find_between(nline, 'step', 'out of').strip(' ')
        t_max = find_between(nline, 'T_max = ', ';').strip(' ')
        q_he4_cu = find_between(nline, 'q_he4_cu = ', ';').strip(' ')
        q_cu_he3 = find_between(nline, 'q_cu_he3 = ', ';').strip(' ')
        if step and t_max:
            i += 1
            step_list.append(int(step))
            t_max_list.append(float(t_max) - t_init)
            q_he4_cu_list.append(abs(float(q_he4_cu) - float(q_deposit))/float(q_deposit))
            q_cu_he3_list.append(abs(float(q_cu_he3) - float(q_deposit))/float(q_deposit))
            if i % 10 == 0:
                plot_t_max.set_data(step_list, t_max_list)
                plot_q_he4_cu.set_data(step_list, q_he4_cu_list)
                plot_q_cu_he3.set_data(step_list, q_cu_he3_list)
                plt.xlim([0, int(step)])
                plt.draw()
                plt.pause(0.00001)
        if 'step' in nline:
            print(Fore.WHITE + nline)
        else:
            if 'relaxation' in nline:
                print(Fore.MAGENTA + nline)
                printted = True
            else:
                if 'updating thermal properties' in nline:
                    print(Fore.YELLOW + nline)
                    printed = True
                else:
                    if 'saving output file' in nline:
                        print(Fore.CYAN + nline)
                        printed = True
                    else:
                        print(Fore.RED + nline)
plt.show()
