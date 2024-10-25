# from fastapi import FastAPI
# from pydantic import BaseModel

# app = FastAPI()

# # Example data model for POST request body
# class ExampleData(BaseModel):
#     name: str
#     value: int

# # Route to handle GET requests at root '/'
# @app.get("/")
# def read_root():
#     return {"message": "Hello World"}  # Or change to {"status": "OK"}

# # Route to handle GET requests
# @app.get("/get_example")
# def get_example():
#     return {"message": "This is a GET request example"}

# # Route to handle POST requests
# @app.post("/post_example")
# def post_example(data: ExampleData):
#     return {"received_data": data}


# from fastapi import FastAPI
# from pydantic import BaseModel
# from typing import Dict

# app = FastAPI()

# # Updated Example data model for POST request body to include JSON data
# class ExampleData(BaseModel):
#     name: str
#     value: int
#     data: Dict  # This field accepts a JSON object

# # Route to handle GET requests at root '/'
# @app.get("/")
# def read_root():
#     return {"message": "Hello World"}

# # Route to handle GET requests
# @app.get("/get_example")
# def get_example():



#     return {"message": "This is a GET request example"}

# # Route to handle POST requests with a JSON object
# @app.post("/post_example")
# def post_example(data: ExampleData):
#     return {
#         "received_name": data.name,
#         "received_value": data.value,
#         "received_json_data": data.data
#     }













# from fastapi import FastAPI, File, UploadFile
# import shutil
# from pathlib import Path
# from Bio.Blast.Applications import NcbiblastxCommandline
# from Bio.Blast import NCBIXML
# import os

# app = FastAPI()

# @app.get("/")
# def read_root():
#     return {"message": "Hello, FastAPI!"}


# @app.post("/run-blast/")
# async def run_blast(file: UploadFile = File(...)):
#     # Save the uploaded file
#     file_location = f"./{file.filename}"
#     with open(file_location, "wb+") as file_object:
#         shutil.copyfileobj(file.file, file_object)

#     # Define filenames for output
#     file_input = file_location
#     file_output = f"{file_location}_blastResults.xml"
#     final_file = f"{file_location}_finalFile.csv"

#     # Run the BLAST search using Biopython's NcbiblastxCommandline
#     blastx_cline = NcbiblastxCommandline(
#         query=file_input, 
#         db="./all/nr_all_amr", 
#         evalue=0.001, 
#         outfmt=5,  # XML output format
#         max_target_seqs=1000,
#         num_threads=5,
#         out=file_output
#     )

#     try:
#         stdout, stderr = blastx_cline()  # Running the command
#     except Exception as e:
#         return {"error": f"Error running BLAST: {str(e)}"}

#     # Parse the BLAST results using Biopython
#     with open(file_output) as result_handle:
#         blast_records = NCBIXML.parse(result_handle)
#         result_data = []
#         for blast_record in blast_records:
#             for alignment in blast_record.alignments:
#                 for hsp in alignment.hsps:
#                     # Extract relevant data from the BLAST output
#                     result_data.append({
#                         "query": blast_record.query,
#                         "subject": alignment.title,
#                         "e_value": hsp.expect,
#                         "identity": hsp.identities,
#                         "length": alignment.length,
#                         "query_seq": hsp.query,
#                         "subject_seq": hsp.sbjct
#                     })

#     # Save the processed results into a CSV file
#     import csv
#     with open(final_file, mode='w', newline='') as csv_file:
#         fieldnames = ['query', 'subject', 'e_value', 'identity', 'length', 'query_seq', 'subject_seq']
#         writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

#         writer.writeheader()
#         for row in result_data:
#             writer.writerow(row)

#     # Optionally delete the intermediate BLAST result file
#     Path(file_output).unlink()

#     # Return the result file path or processed result
#     return {"message": f"BLAST search completed. Final results saved to {final_file}"}

# To run the server:
# uvicorn your_script_name:app --reload






























from fastapi import FastAPI, File, UploadFile
import shutil
from pathlib import Path
import os
from vasoula_1component import find_amr_genes_module
import pandas as pd
from fastapi.middleware.cors import CORSMiddleware


app = FastAPI()

# Add CORS middleware to allow requests from your Next.js frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # Allows Next.js running on port 3000
    allow_credentials=True,
    allow_methods=["*"],  # Allow all HTTP methods (POST, GET, etc.)
    allow_headers=["*"],  # Allow all headers
)

@app.get("/")
def read_root():
    return {"message": "Hello, FastAPI!"}

@app.post("/run-blast/")
async def run_blast(
    file: UploadFile = File(...),
    db_path: str = "vasoula_1component/databases/all/all_amr",
    num_threads: int = 4,
    protein_annotation_file: str = "vasoula_1component/databases/genes_annotation_databases.csv",
    protein_start_filter: int = 20,
    identity_threshold: float = 70.0,
    coverage_threshold: float = 70.0,
    output_csv_file: str = "final_results.csv"
):
    # Save the uploaded file (query file)
    file_location = f"./{file.filename}"
    with open(file_location, "wb+") as file_object:
        shutil.copyfileobj(file.file, file_object)

    query_file = file_location  # The query file path
    
    # Run BLASTX and store the results in a DataFrame using the provided module
    try:
        blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(
            query_file, db_path, num_threads
        )
    except Exception as e:
        return {"error": f"Error running BLAST: {str(e)}"}

    # Process the BLAST results using the provided module
    try:
        final_results = find_amr_genes_module.process_blast_results(
            protein_annotation_file, blast_results_df, query_file, 
            protein_start_filter, identity_threshold, coverage_threshold
        )
    except Exception as e:
        return {"error": f"Error processing BLAST results: {str(e)}"}

    # Save the final results to a CSV file
    try:
        final_results.to_csv(output_csv_file, index=False)
    except Exception as e:
        return {"error": f"Error saving final results: {str(e)}"}

    # Optionally, delete the query file after processing
    Path(query_file).unlink()

    return {"message": f"BLAST search and processing completed. Final results saved to {output_csv_file}"}


