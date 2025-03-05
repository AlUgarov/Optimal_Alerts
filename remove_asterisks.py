import os
import re

# Define the main LaTeX file
main_tex_file = "main.tex"  # Change if your main file has a different name

# Patterns to match content to be modified
input_pattern = re.compile(r'\\input\{([^}]+)\}')
asterisk_pattern = re.compile(r'\*+')
footnote_pattern = re.compile(
    r'\\footnotesize\s*\\sym\{\*\}.*?\\sym\{\*\*\}.*?\\sym\{\*\*\*\}.*?\}',
    re.DOTALL
)

# Function to process LaTeX content (removes asterisks and footnotes)
def clean_latex_content(content):
    content = asterisk_pattern.sub("", content)
    content = footnote_pattern.sub("", content)
    return content

# Function to process and save an external table file
def process_table_file(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    cleaned_content = clean_latex_content(content)

    new_file_path = file_path.replace(".tex", "_noast.tex")
    with open(new_file_path, "w", encoding="utf-8") as f:
        f.write(cleaned_content)
    
    return new_file_path

# Process the main LaTeX document
with open(main_tex_file, "r", encoding="utf-8") as f:
    main_content = f.read()

# Apply modifications to the main LaTeX content
cleaned_main_content = clean_latex_content(main_content)

# Find and process all tables included via \input{}
matches = input_pattern.findall(main_content)
updated_files = {}

for table_file in matches:
    table_path = table_file if table_file.endswith(".tex") else table_file + ".tex"

    if os.path.exists(table_path):
        new_table_path = process_table_file(table_path)
        updated_files[table_file] = new_table_path.replace(".tex", "")

# Update \input{} references in the cleaned main content
for old, new in updated_files.items():
    cleaned_main_content = cleaned_main_content.replace(f"\\input{{{old}}}", f"\\input{{{new}}}")

# Save the modified main LaTeX document
new_main_tex_file = main_tex_file.replace(".tex", "_noast.tex")
with open(new_main_tex_file, "w", encoding="utf-8") as f:
    f.write(cleaned_main_content)

print(f"Updated LaTeX document saved as {new_main_tex_file}")