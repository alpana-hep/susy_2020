import os

def search_string_in_files(directory, search_string):
    found_files_with_string = []
    found_files_without_string = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    all_lines_contain_string = True
                    for line_number, line in enumerate(f, 1):
                        if search_string not in line:
                            all_lines_contain_string = False
                            break
                    if all_lines_contain_string:
                        found_files_with_string.append(file_path)
                    else:
                        found_files_without_string.append(file_path)
                        # Replace '/store/user' with 'root://cmseos.fnal.gov//store/user'
                        with open(file_path, 'r') as f:
                            file_content = f.read()
                            modified_content = file_content.replace('/store/user', 'root://cmseos.fnal.gov//store/user')
                        with open(file_path, 'w') as f:
                            f.write(modified_content)
    return found_files_with_string, found_files_without_string

directory_to_search = './'  # Specify the directory you want to search in
search_string = 'root://cmseos.fnal.gov/'

found_files_with_string, found_files_without_string = search_string_in_files(directory_to_search, search_string)

print("Files where all lines contain the search string:")
if found_files_with_string:
    for file_path in found_files_with_string:
        print(file_path)
else:
    print("No files were found where all lines contain the search string.")

print("\nFiles where not all lines contain the search string:")
if found_files_without_string:
    for file_path in found_files_without_string:
        print(file_path)
else:
    print("All files contain the search string in at least one line.")
