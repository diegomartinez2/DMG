# Import the xlrd module
import xlrd

# Open the Workbook
workbook = xlrd.open_workbook('sales.xlsx')

# Open the worksheet
worksheet = workbook.sheet_by_index(0)

# Iterate the rows and columns
for i in range(0, 5):
    for j in range(0, 3):
# Print the cell values with tab space
        print(worksheet.cell_value(i, j), end='\t')
    print('')

#------------with pandas---------------
# Import pandas
import pandas as pd

# Load the xlsx file
excel_data = pd.read_excel('sales.xlsx')
# Read the values of the file in the dataframe
data = pd.DataFrame(excel_data, columns=[
‘Sales Date’, ‘Sales Person’, ‘Amount’])
# Print the content
print(“The content of the file is:\n”, data)
