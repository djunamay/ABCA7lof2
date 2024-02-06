# list all files in synapse submission for manifest

find "./synapse" -type f | while read -r file; do
    echo "$file" >> "all_files.txt"
done
