import os
import glob
import xlrd

# Get the current directory
current_directory = os.getcwd()

# Find all .xlsx files in the current directory
excel_files = glob.glob(os.path.join(current_directory, "*.xlsx"))

# Iterate over each Excel file
for file_name in excel_files:
    # Open the workbook
    workbook = xlrd.open_workbook(file_name)

    # Open the worksheet (assuming the first sheet is the one to be processed)
    worksheet = workbook.sheet_by_index(0)

    # Iterate the rows and columns
    for i in range(worksheet.nrows):
        for j in range(worksheet.ncols):
            # Print the cell values with tab space
            print(worksheet.cell_value(i, j), end='\t')
        print('')
