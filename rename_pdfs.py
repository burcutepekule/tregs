#!/usr/bin/env python3
"""
Script to rename PDF files with their article titles.
Extracts the title from the PDF metadata or content and renames the file.
"""

import os
import re
from pathlib import Path
import PyPDF2

def clean_filename(title):
    """
    Clean the title to make it a valid filename.
    Removes special characters and limits length.
    """
    # Remove or replace invalid filename characters
    title = re.sub(r'[<>:"/\\|?*]', '', title)
    # Replace multiple spaces with single space
    title = re.sub(r'\s+', ' ', title)
    # Remove leading/trailing whitespace
    title = title.strip()
    # Limit length to 100 characters
    if len(title) > 100:
        title = title[:100].rsplit(' ', 1)[0]  # Cut at last space before 100 chars
    return title

def extract_title_from_pdf(pdf_path):
    """
    Extract title from PDF metadata or first page content.
    """
    try:
        with open(pdf_path, 'rb') as file:
            pdf_reader = PyPDF2.PdfReader(file)
            
            # First, try to get title from metadata
            metadata = pdf_reader.metadata
            if metadata and metadata.title:
                title = metadata.title.strip()
                if title and len(title) > 5:  # Reasonable title length
                    return title
            
            # If no metadata title, extract from first page
            if len(pdf_reader.pages) > 0:
                first_page = pdf_reader.pages[0]
                text = first_page.extract_text()
                
                # Get first non-empty line as potential title
                lines = [line.strip() for line in text.split('\n') if line.strip()]
                if lines:
                    # Often the title is in the first few lines
                    # Filter out very short lines (likely not titles)
                    for line in lines[:10]:
                        if len(line) > 15 and not line.isupper():
                            return line
                    # Fallback to first line if nothing better found
                    return lines[0]
    
    except Exception as e:
        print(f"Error reading {pdf_path.name}: {e}")
        return None
    
    return None

def rename_pdfs_in_folder(folder_path, dry_run=True):
    """
    Rename all PDF files in the specified folder.
    
    Args:
        folder_path: Path to the folder containing PDFs
        dry_run: If True, only show what would be renamed without actually renaming
    """
    folder = Path(folder_path)
    
    if not folder.exists():
        print(f"Error: Folder {folder_path} does not exist")
        return
    
    pdf_files = list(folder.glob("*.pdf"))
    
    if not pdf_files:
        print("No PDF files found in the folder")
        return
    
    print(f"Found {len(pdf_files)} PDF files")
    print("=" * 80)
    
    renamed_count = 0
    failed_count = 0
    
    for pdf_file in sorted(pdf_files):
        print(f"\nProcessing: {pdf_file.name}")
        
        title = extract_title_from_pdf(pdf_file)
        
        if not title:
            print(f"  ⚠️  Could not extract title")
            failed_count += 1
            continue
        
        clean_title = clean_filename(title)
        new_filename = f"{clean_title}.pdf"
        new_path = pdf_file.parent / new_filename
        
        # Check if new filename already exists (and is different from current)
        if new_path.exists() and new_path != pdf_file:
            print(f"  ⚠️  Target filename already exists: {new_filename}")
            failed_count += 1
            continue
        
        # Skip if the filename is essentially the same
        if pdf_file.name == new_filename:
            print(f"  ℹ️  Already has correct name")
            continue
        
        print(f"  📄 Title: {title[:80]}")
        print(f"  ➜  New name: {new_filename}")
        
        if not dry_run:
            try:
                pdf_file.rename(new_path)
                print(f"  ✓  Renamed successfully")
                renamed_count += 1
            except Exception as e:
                print(f"  ✗  Error renaming: {e}")
                failed_count += 1
        else:
            renamed_count += 1
    
    print("\n" + "=" * 80)
    print(f"Summary:")
    print(f"  - Would rename: {renamed_count} files" if dry_run else f"  - Renamed: {renamed_count} files")
    print(f"  - Failed/Skipped: {failed_count} files")
    
    if dry_run:
        print("\nThis was a DRY RUN. No files were actually renamed.")
        print("Run with dry_run=False to perform actual renaming.")

if __name__ == "__main__":
    import sys
    
    # Get folder path from command line argument or use default
    if len(sys.argv) > 1:
        folder_path = sys.argv[1]
    else:
        # Default path based on the screenshot
        folder_path = "~/papers"
    
    # Expand user home directory
    folder_path = os.path.expanduser(folder_path)
    
    # Check if user wants to actually rename (use --rename flag)
    dry_run = "--rename" not in sys.argv
    
    print("PDF Renaming Script")
    print("=" * 80)
    print(f"Target folder: {folder_path}")
    print(f"Mode: {'DRY RUN (no changes will be made)' if dry_run else 'ACTUAL RENAMING'}")
    print()
    
    rename_pdfs_in_folder(folder_path, dry_run=dry_run)
