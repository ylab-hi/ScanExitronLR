#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 20:44:36 2022

@author: jpfry
"""
__version__ = 'v.2'
import os, sys, subprocess

def main():

    file_abs_path = os.path.abspath(os.path.dirname(__file__))
    crtAbsPath = os.getcwd()

    try:
        task = sys.argv[1]
        if task not in ['extract', 'annotate']:
            if task not in ['-v', '--version']: print(f'\nERROR: Unknown usage.')
            raise ValueError
        else:
            if task == 'extract':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(['python', 'extract.py'] + args, check = True)
            elif task == 'annotate':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(['python', 'annotate.py'] + args, check = True)
    except:
        print(f'\nProgram:\tScanExitronLR')
        print(f'Version:\t{__version__}')
        print(f'Usage:\t\tselr <command> [options]')
        print(f'\nCommands:\textract\t\t\tTool for extracting exitrons and creating exitron files')
        print(f'\t\tannotate\t\tTool for annotating exitrons from exitron files')


if __name__ == '__main__':
    main()
