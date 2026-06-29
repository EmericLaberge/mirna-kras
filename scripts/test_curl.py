"""One-shot smoke test of the RIMap-RISC API.

Sends a single POST request for the SIRT1-201 / miR-34a-5p pair and prints the
JSON response. Use it to verify the API is reachable and your network is OK
before launching the full ``rimap_pipeline.py`` run.
"""
import requests
import json


# URL of the RIMap-RISC API datapoint
url = 'https://rimap-risc.api.major.iric.ca/api/data'

# JSON payload for the request
payload = {
   "gencode": "v46",
   "transcriptName": "SIRT1-201",
   "miRName": "miR-34a-5p",
   "minSeed": "6mer",
   "UTR5": False,
   "CDS": True,
   "UTR3": True,
   "supplementaryPaired": False,
   "seedAccessible": 0.0,
   "supplementaryAccessible": 0.0,
   "conserved": 0
}

try:
   # Send the POST request
   response = requests.post( url, json=payload )
   # Check if the request was successful
   if response.status_code == 200:
       print( "Request was successful!" )
       response_data = response.json()
       # Print the entire response as a formatted JSON string
       print( json.dumps( response_data, indent=4 ) )
       # Extract specific parts of the response
       miRNA_target_pair = response_data.get( 'miRNATargetPair', 'N/A' )
       binding_landscape = response_data.get( 'bindingLandscape', 'No interactions found' )
       transcript_regions = response_data.get( 'transcriptRegions', [] )
       binding_positions = response_data.get( 'bindingPositions', [] )
       transcript_accessibility = response_data.get( 'transcriptAccessibility', [] )
       print( "\n--- Interaction Results ---" )
       print( f"miRNA-Target Pair: {miRNA_target_pair}" )
       print( f"Binding Landscape:\n{binding_landscape}" )
       print( f"Transcript Regions: {transcript_regions}" )
       print( f"Binding Positions: {binding_positions}" )
       print( f"Accessibility Scores: {transcript_accessibility}" )
   else:
       print( f"Request failed with status code: {response.status_code}" )
       print( response.text )
except Exception as e:
   print( f"An error occurred: {e}" )

# CURL:
# curl -X POST "https://rimap-risc.api.major.iric.ca/api/data" \
#   -H "Content-Type: application/json" \
#   -d '{"gencode": "v46", "transcriptName": "SIRT1-201", "miRName": "miR-34a-5p", "UTR3": true}'
