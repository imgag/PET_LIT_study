import os
import sys
import shlex
import pandas as pd
import subprocess


def vcf_add_filter_list(filter_column: str, add_filter_list: list):

    for f in add_filter_list:
        filter_column = vcf_add_filter(filter_column, f)

    return filter_column

def vcf_add_filter(filter_column: str, add_filter: str):

    if filter_column == 'PASS':
        filter_column = ''

    if pd.isnull(filter_column):
        filter_column = ''

    splitted = filter_column.split(';')
    splitted.append(add_filter)

    # remove empty elements
    splitted = list(filter(None, splitted))

    return ';'.join(splitted)

def run_python_tool(command):

    cmd = '{} {}'.format(sys.executable,command)
    print(cmd)

    try:
        if sys.platform=='linux':
            cmd = shlex.split(cmd)
            output = subprocess.run(cmd, capture_output=True, check=True)
        elif sys.platform=='win32':
            CREATE_NO_WINDOW = 0x08000000
            output = subprocess.run(cmd, capture_output=True, check=True, creationflags=CREATE_NO_WINDOW)
    except subprocess.CalledProcessError as exc:
        print(cmd)
        print(subprocess.CalledProcessError)
        print("Status : FAIL", exc.returncode, exc.stderr)
        sys.exit(0)

    return output