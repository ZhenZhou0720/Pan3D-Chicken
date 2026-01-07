import sys
import numpy as np
import pandas as pd
import bbi

inbw = sys.argv[1]
region_file = sys.argv[2]
nbins = int(sys.argv[3])
output_file = sys.argv[4] if len(sys.argv) >=5 else "average_signal.csv"

df = pd.read_csv(region_file)
stackup = bbi.stackup(inbw, df['chrom'], df['start'], df['end'], bins=nbins)
average_signal = np.nanmean(stackup, axis=0)

result_df = pd.DataFrame({
    "bin": range(1, nbins+1),
    "average_signal": average_signal
})
result_df.to_csv(output_file, index=False)