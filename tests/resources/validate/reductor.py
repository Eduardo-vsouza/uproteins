import pandas as pd
import sys


def filter_scans(results_04, pin_file, output):
    results = pd.read_csv(results_04, sep='\t')
    scans  = results["ScanNum"].tolist()
    denovo = results["DeNovoScore"].tolist()
    pin = pd.read_csv(pin_file, sep='\t')
    pin = pin[pin["ScanNr"].isin(scans)]
    pin = pin[pin["DeNovoScore"].isin(denovo)]
    pin.to_csv(output, sep='\t', index=False)
    

if __name__ == '__main__':
    filter_scans(sys.argv[1], sys.argv[2], sys.argv[3])
