#!/usr/bin/env python3
"""
Created on Wed Mar 16 20:44:36 2022

@author: jpfry
"""
__version__ = 'v1.1.6'
import os
import sys
import subprocess


def main():

    file_abs_path = os.path.abspath(os.path.dirname(__file__))

    try:
        task = sys.argv[1]
        if task not in ['extract', 'annotate']:
            if task not in ['-v', '--version']:
                print(f'\nERROR: Unknown usage.')
            raise ValueError
        else:
            if task == 'extract':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(
                    ['python', f'{file_abs_path}/extract.py'] + args, check=True)
            elif task == 'annotate':
                print('\n')
                args = sys.argv[2:]
                subprocess.run(
                    ['python', f'{file_abs_path}/annotate.py'] + args, check=True)
    except:
        print(f'\nProgram:\tScanExitronLR')
        print(f'Version:\t{__version__}')
        print(f'Usage:\t\tselr <command> [options]')
        print(f'\nCommands:\textract\t\t\tTool for extracting exitrons and creating exitron files')
        print(f'\t\tannotate\t\tTool for annotating exitrons from exitron files')


if __name__ == '__main__':
    main()
