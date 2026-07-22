$baseUrl = "http://localhost:1234/v1"

# Load new filtered network
Invoke-RestMethod -Uri "$baseUrl/commands/network/load file" -Method Post -Body '{"file": "d:/Zika_wetlab/meta-analysis/results/STRING_PPI_Network_Cytoscape.graphml"}' -ContentType 'application/json' | Out-Null
Start-Sleep -Seconds 2

# Get the latest network ID
$networkId = (Invoke-RestMethod -Uri "$baseUrl/networks")[0]

# Apply the flat reference style
Invoke-RestMethod -Uri "$baseUrl/apply/styles/ReferenceImageStyle/$networkId" -Method Get | Out-Null

# Apply unclustered force-directed layout
Invoke-RestMethod -Uri "$baseUrl/commands/layout/force-directed" -Method Post -Body '{}' -ContentType "application/json" | Out-Null
Start-Sleep -Seconds 3

# Export final image
Invoke-RestMethod -Uri "$baseUrl/networks/$networkId/views/first.png?h=2500" -OutFile "d:\Zika_wetlab\meta-analysis\plots\12_Cytoscape_Reference_Filtered.png"
