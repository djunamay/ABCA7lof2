# copy lipidomics data

for file in "../raw_data/mass_spec_core_raw_data/SUB12418"/*.raw; do
    echo "$file"
    filename=$(basename -- "$file")
    cp "$file" "./synapse/lipidomics/$filename"
done
