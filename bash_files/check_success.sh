
for file in ../raw_data/cellranger_counts_out/*.out
do
    if grep -q 'Pipestance completed successfully!' "$file"
    then
      echo $file + 'success'
    else
      echo $file + 'FAILED'
    fi
done
