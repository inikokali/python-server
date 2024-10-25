import filecmp

file1 = "final_results.csv"
file2 = "vasoula_1component/results.csv"

# Compare the two files
are_files_equal = filecmp.cmp(file1, file2, shallow=False)

if are_files_equal:
    print("The files are the same.")
else:
    print("The files are different.")
