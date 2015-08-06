#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# created: 05.06.2015
# author: Aleksey Komissarov
# contact: ad3002@gmail.com

import argparse
import time
import sys
import subprocess


def run_asap(cmd_list, cpu_num=10, mock=False):
    """
    Run large number of cmd_list in parallel with subprocess.Popen
    :param cmd_list: list of cmd_list for shell
    :param cpu_num: number of cpu_num
    :param mock: don't run command
    :type cmd_list: list
    :type cpu_num: int
    :type mock: bool
    """
    if not isinstance(cmd_list, list):
        message = "Expected list get %s" % str(cmd_list)
        print(message)
        raise Exception(message)
    running = []
    while cmd_list:
        while len(running) > cpu_num:
            print("Checking %s processes from (%s)" % (len(running), len(cmd_list)))
            for j, p in enumerate(running):
                returncode = p.poll()
                if returncode is not None:
                    if returncode == 0:
                        print('A process returned: %s (remains %s)' %
                              (p.returncode, len(cmd_list)))
                    else:
                        print('A process returned error: %s (remains '
                              '%s)' % (p.returncode, len(cmd_list)))
                    running[j] = None
                    running = [x for x in running if x is not None]
                    break
            time.sleep(1)
        curr_cmd = cmd_list.pop()
        print(curr_cmd)
        if not mock:
            running.append(subprocess.Popen(curr_cmd, shell=True))

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Run Cookiecutter for multiple input files in '
                    'parallel.')
    parser.add_argument('-i', '--input', help='Comma-separated list '
                                              'of single-end read '
                                              'FASTQ files')
    parser.add_argument('-1', '--fastq1',
                        help='Comma-separated list of left FASTQ files')
    parser.add_argument('-2', '--fastq2', 
                        help='Comma-separated list of right FASTQ '
                             'files')
    parser.add_argument('-s', '--single_end', action='store_true',
                        help='input data are single-end reads',
                        required=True)
    parser.add_argument('-c', '--command', 
                        help='Cookiecutter subprogram', required=True)
    parser.add_argument('-o', '--out', help='Output folder', 
                        required=True)
    parser.add_argument('-f', '--fragments', help='Kmer library', 
                        required=True)
    parser.add_argument('-P', '--cpus', help='CPUs', required=False, 
                        default=4)
    parser.add_argument('-g', '--polyG',
                        help='Length of polyG/polyC track to filter '
                             'out', required=False, default=23)
    parser.add_argument('-l', '--length', help='Minimal read length',
                        required=False, default=50)
    parser.add_argument('-d', '--dustcutoff', 
                        help='Cutoff for DUST algorithm',
                        required=False, default=None)
    parser.add_argument('-k', '--dustk', help='K for DUST algorithm',
                        required=False, default=None)
    parser.add_argument('-q', '--mq', help='Mean read quality',
                        required=False, default=20)
    args = vars(parser.parse_args())

    fastq_files1 = fastq_files2 = fastq_files = []

    if args['single_end']:
        fastq_files = args["input"].split(",")
    else:
        fastq_files1 = args["fastq1"].split(",")
        fastq_files2 = args["fastq2"].split(",")

    cpu = int(args["cpus"])
    out_folder = args["out"]
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

    available_commands = ('remove', 'rm_reads', 'separate', 'extract')
    if command not in available_commands:
        print("Unknown command. Available: %s" %
              ', '.join(available_commands))
        sys.exit(2)

    options = {
        "command": command,
        "out": out_folder,
        "fragments": fragments,
        }
    command_options = "%(command)s -o %(out)s --fragments %(" \
                      "fragments)s" % options
    commands = []
    for i in range(len(fastq_files1)):
        if args["single_end"]:
            command_data = "-i {}".format(fastq_files[i])
        else:
            command_data = "-1 {} -2 {}".format(fastq_files1[i],
                                                fastq_files2[i])
        if command == "rm_reads":
            command_options += '--polyG {}'.format(polyG)
            command_options += '--length {}'.format(length)
            command_options += '--mq {}'.format(mq)
            if dustk and dustcutoff:
                command_options += '--dust_cutoff {}'.format(dustcutoff)
                command_options += '--dust_k {}'.format(dustk)
        commands.append(command_options + command_data)

    run_asap(commands, cpu_num=10, mock=False)
