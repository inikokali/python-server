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





from fastapi import FastAPI, File, UploadFile, Form, HTTPException
import shutil
from pathlib import Path
import os
from vasoula_1component import find_amr_genes_module
import pandas as pd
from fastapi.middleware.cors import CORSMiddleware
from Bio.Blast.Applications import NcbiblastxCommandline
import requests
import json  # Add this line
from fastapi.responses import FileResponse
import os


app = FastAPI()

# Add CORS middleware to allow requests from your Next.js frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Hello, FastAPI!"}


@app.post("/run-blast/")
async def run_blast(
    file: UploadFile = File(...),
    jobId: str = Form(...),
    db_path: str = "vasoula_1component/databases/all/all_amr",
    num_threads: int = 4,
    protein_annotation_file: str = "vasoula_1component/databases/genes_annotation_databases.csv",
    protein_start_filter: int = 20,
    identity_threshold: float = 70.0,
    coverage_threshold: float = 70.0,
):
    try:
        # Save the uploaded file locally
        file_location = f"./{file.filename}"
        with open(file_location, "wb+") as file_object:
            shutil.copyfileobj(file.file, file_object)

        # Check if file is empty
        if os.path.getsize(file_location) == 0:
            raise HTTPException(status_code=400, detail="Uploaded file is empty")

        # Additional check to ensure file contains valid content for BLASTX
        with open(file_location, "r") as f:
            lines = f.readlines()
            if not any(line.startswith(">") for line in lines):
                raise HTTPException(status_code=400, detail="Uploaded file contains no valid FASTA sequences")

        # Define output file path for the BLAST results
        output_csv_file = f"/home/vkostoula/argis/python-server/final_results_{jobId}.csv"
        temp_blast_output = f"{file_location}_blastResults.xml"  # Temporary file for BLAST output

        # Run the BLASTX search using Biopython's NcbiblastxCommandline
        try:
            blastx_cline = NcbiblastxCommandline(
                query=file_location,
                db=db_path,
                evalue=0.001,
                outfmt=5,  # XML output format
                max_target_seqs=1000,
                num_threads=num_threads,
                out=temp_blast_output
            )
            stdout, stderr = blastx_cline()  # Run the BLASTX command
            print("BLASTX process completed successfully")

        except Exception as e:
            print("Error in BLASTX process:", e)
            raise HTTPException(status_code=500, detail=f"Error running BLASTX: {str(e)}")

        # Process the BLAST results with the find_amr_genes_module
        try:
            # Read BLAST output and store results in a DataFrame
            blast_results_df = find_amr_genes_module.run_blastx_to_dataframe(
                file_location, db_path, num_threads
            )

            # Further process the DataFrame to filter and annotate results
            final_results = find_amr_genes_module.process_blast_results(
                protein_annotation_file, blast_results_df, file_location, 
                protein_start_filter, identity_threshold, coverage_threshold
            )

            # Save the final processed results to a CSV file
            final_results.to_csv(output_csv_file, index=False)
            print(f"Final results saved to {output_csv_file}")

        except Exception as e:
            print("Error in processing BLAST results:", e)
            raise HTTPException(status_code=500, detail=f"Error processing BLAST results: {str(e)}")
        
        
        try:
            update_url = "http://localhost:3000/api/jobStatus"
            response = requests.post(update_url, json={
                "jobId": jobId,
                "status": "completed",
                "resultPath": output_csv_file
            })
            response.raise_for_status()  # Raise error if the update fails
            print(f"Job status updated to 'completed' for jobId: {jobId}")
        except Exception as e:
            print(f"Error updating job status to completed: {e}")

        # Clean up temporary files
        Path(temp_blast_output).unlink()
        Path(file_location).unlink()

        return {"message": f"Job {jobId} completed. Results saved to {output_csv_file}"}

    except Exception as e:
        print("General Error:", e)
        raise HTTPException(status_code=500, detail="An unexpected error occurred")



# New endpoint to serve CSV files based on job ID
@app.get("/get-result/{job_id}")
async def get_result(job_id: str):
    file_path = f"/home/vkostoula/argis/python-server/final_results_{job_id}.csv"
    
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Result file not found")
    
    return FileResponse(path=file_path, media_type="text/csv", filename=f"final_results_{job_id}.csv")

