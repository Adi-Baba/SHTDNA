
import os
import glob
import sys

# Try to import pypdf or PyPDF2
try:
    import pypdf
    print("Using pypdf")
except ImportError:
    try:
        import PyPDF2 as pypdf
        print("Using PyPDF2")
    except ImportError:
        print("No PDF library found. Attempting to install pypdf...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pypdf"])
        import pypdf
        print("Installed and using pypdf")

directory = r"d:\OnlyHST\SHTDNA\SHT_DNA"
pdf_files = glob.glob(os.path.join(directory, "*.pdf"))

if not pdf_files:
    print(f"No PDF files found in {directory}")
    exit()

for pdf_path in pdf_files:
    print(f"Processing {pdf_path}...")
    try:
        reader = pypdf.PdfReader(pdf_path)
        text = ""
        for page in reader.pages:
            extract = page.extract_text()
            if extract:
                text += extract + "\n"
        
        output_filename = os.path.basename(pdf_path) + ".txt"
        output_path = os.path.join(directory, output_filename)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"Saved to {output_path}")
    except Exception as e:
        print(f"Failed to extract {pdf_path}: {e}")
