#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File_organizer.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
# File organizer, not tested yet by me.
import os
import shutil

# Define source folder to analyze
source_folder = r"C:\Users\USERNAME\Desktop\to_organize"  # Update with your source folder path

# Define destination folder paths
destination_root = r"C:\Users\USERNAME\Desktop\Organized"  # Root for organized files
DESTINATIONS = {
    "Images": os.path.join(destination_root, "Images"),
    "Documents": os.path.join(destination_root, "Documents"),
    "Spreadsheets": os.path.join(destination_root, "Spreadsheets"),
    "Audio": os.path.join(destination_root, "Audio"),
    "Videos": os.path.join(destination_root, "Videos"),
    "Shortcuts": os.path.join(destination_root, "Shortcuts"),
    "Books": os.path.join(destination_root, "Books"),
    "Scripts": os.path.join(destination_root, "Scripts"),
    "Miscellaneous": os.path.join(destination_root, "Miscellaneous"),
    "Zipped": os.path.join(destination_root, "Zipped"),
    "Folders": os.path.join(destination_root, "Folders"),
}

# File type extensions for categories
FILE_CATEGORIES = {
    "Images": [".jpg", ".jpeg", ".png", ".gif", ".jfif", ".bmp", ".tiff", ".svg", ".html", ".webp"],
    "Documents": [".doc", ".docx", ".pdf", ".txt", ".ppt", ".pptx"],
    "Spreadsheets": [".xls", ".xlsx", ".csv", ".json", ".xml"],
    "Audio": [".mp3", ".wav", ".aac", ".flac", ".ogg", ".wma"],
    "Videos": [".mp4", ".mkv", ".avi", ".mov", ".wmv", ".flv"],
    "Scripts": [".ps1", ".py"],
    "Shortcuts": [".lnk", ".url", ".webloc"],
    "Zipped": [".zip", ".rar", ".7z"],
    "Books": [".epub"]
}

# Create destination folders if they don't exist
for folder in DESTINATIONS.values():
    os.makedirs(folder, exist_ok=True)

def organize_files(source):
    for item in os.listdir(source):
        item_path = os.path.join(source, item)

        # If the item is a directory, move it to "Folders"
        if os.path.isdir(item_path):
            shutil.move(item_path, os.path.join(DESTINATIONS["Folders"], item))
            print(f"Moved folder: {item} to Folders")
            continue

        # If the item is a file, process it
        if os.path.isfile(item_path):
            _, extension = os.path.splitext(item)
            extension = extension.lower()

            # Check the extension against file categories
            for category, extensions in FILE_CATEGORIES.items():
                if extension in extensions:
                    shutil.move(item_path, os.path.join(DESTINATIONS[category], item))
                    print(f"Moved file: {item} to {category}")
                    break
            else:
                # If no category matches, move to "Miscellaneous"
                shutil.move(item_path, os.path.join(DESTINATIONS["Miscellaneous"], item))
                print(f"Moved file: {item} to Miscellaneous")

if __name__ == "__main__":
    if not os.path.exists(source_folder):
        print(f"The source folder '{source_folder}' does not exist. Please check the path.")
    else:
        organize_files(source_folder)
        print("Organization complete!")
