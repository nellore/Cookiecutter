#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2015
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import argparse
import os
import sys


def run_asap(commands, cpu=10, mock=False):
    """
    Run large number of commands in parallel with subprocess.Popen
    :param commands: list of commands for shell
    :param cpu: number of cpu
    :param mock: don't run command
    """
    if not isinstance(commands, list):
        message = "Expected list get %s" % str(commands)
        print(message)
        raise Exception(message)
    running = []
    while commands:
        while len(running) > cpu:
            print("Checking %s processes from (%s)" % (len(running), len(commands)))
            for i, p in enumerate(running):
                returncode = p.poll()
                if returncode is not None:
                    if returncode == 0:
                        print('A process returned: %s (remains %s)' % (p.returncode, len(commands)))
                    else:
                        print('A process returned error: %s (remains %s)' % (p.returncode, len(commands)))
                    running[i] = None
                    running = [x for x in running if x is not None]
                    break
            time.sleep(1)
        command = commands.pop()
        print(command)
        if not mock:
            running.append(subprocess.Popen(command, shell=True))

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Run analysis in parallel.')
    parser.add_argument('-1','--fastq1', help='Comma-separated list of left fastq files', required=True)
    parser.add_argument('-2','--fastq2', help='Comma-separated list of right fastq files', required=True)
    parser.add_argument('-c','--command', help='Cookiecutter subprogram', required=True)
    parser.add_argument('-o','--out', help='Output folder', required=True)
    parser.add_argument('-f','--fragments', help='Kmer library', required=True)
    parser.add_argument('-P','--cpus', help='CPUs', required=False, default=4)
    parser.add_argument('-g','--polyG', help='Length of polyG/polyC track to filter out', required=False, default=23)
    parser.add_argument('-l','--legnth', help='Minimal read length', required=False, default=50)
    parser.add_argument('-d','--dustcutoff', help='Cutoff for DUST algorithm', required=False, default=None)
    parser.add_argument('-k','--dustk', help='K for DUST algorithm', required=False, default=None)
    parser.add_argument('-q','--mq', help='Mean read quality', required=False, default=20)
    args = vars(parser.parse_args())

    fastq_files1 = args["fastq1"].split(",")
    fastq_files2 = args["fastq2"].split(",")
    cpu = int(args["cpus"])
    out_folder = int(args["out"])
    command = args["command"]
    fragments = args["fragments"]
    polyG = args["polyG"]
    length = args["length"]
    dustk = args["dustk"]
    dustcutoff = args["dustcutoff"]
    mq = args["mq"]

    if not len(fastq_files1) == len(fastq_files2):
        print("Not equal number of left and right fastq files")
        sys.exit(2)

    available_commands = "remove, rm_reads, separate, extract"
    if not command in available_commands:
        print("Unknow command. Avaiable: %s" % available_commands)
        sys.exit(2)

    commands = []
    for i in range(len(fastq_files1)):
        data = {
            "command": command,
            "left": fastq_files1[i],
            "right": fastq_files2[i],
            "out": out_folder,
            "fragments": fragments,
            "polyG": polyG,
            "length": length,
            "dustcutoff": dustcutoff,
            "dustk": dustk,
            "mq": mq,
        }
        if command == "rm_reads":
            if dustk and dustcutoff:
                command = "%(command)s -1 %(left)s -2 %(right)s -o %(out)s --fragments %(fragments)s  --polyG %(polyG)s --length %(length)s --dust_cutoff %(dustcutoff)s --dust_k %(dustk)s --mq %(mq)s" % data
            else:
                command = "%(command)s -1 %(left)s -2 %(right)s -o %(out)s --fragments %(fragments)s  --polyG %(polyG)s --length %(length)s --mq %(mq)s" % data
        else:   
            command = "%(command)s -1 %(left)s -2 %(right)s -o %(out)s --fragments %(fragments)s" % data
        commands.append(command)

    run_asap(commands, cpu=10, mock=False)

    