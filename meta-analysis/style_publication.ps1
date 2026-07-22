$baseUrl = "http://localhost:1234/v1"

# Create exact match style from reference image
$style = @{
    title = "ReferenceImageStyle"
    defaults = @(
        @{ visualProperty = "NODE_SHAPE"; value = "ROUND_RECTANGLE" },
        @{ visualProperty = "NODE_WIDTH"; value = 65.0 },
        @{ visualProperty = "NODE_HEIGHT"; value = 35.0 },
        @{ visualProperty = "NODE_FILL_COLOR"; value = "#D35252" },
        @{ visualProperty = "NODE_BORDER_WIDTH"; value = 1.0 },
        @{ visualProperty = "NODE_BORDER_PAINT"; value = "#FFFFFF" },
        @{ visualProperty = "NODE_LABEL_COLOR"; value = "#FFFFFF" },
        @{ visualProperty = "NODE_LABEL_FONT_SIZE"; value = 12 },
        @{ visualProperty = "NODE_LABEL_POSITION"; value = "C,C,c,0.00,0.00" },
        @{ visualProperty = "NETWORK_BACKGROUND_PAINT"; value = "#FFFFFF" },
        @{ visualProperty = "EDGE_WIDTH"; value = 4.0 },
        @{ visualProperty = "EDGE_STROKE_UNSELECTED_PAINT"; value = "#BDBDBD" },
        @{ visualProperty = "EDGE_TRANSPARENCY"; value = 255 }
    )
    mappings = @(
        @{
            mappingType = "passthrough"
            mappingColumn = "label"
            mappingColumnType = "String"
            visualProperty = "NODE_LABEL"
        }
    )
}

$body = $style | ConvertTo-Json -Depth 10
Invoke-RestMethod -Uri "$baseUrl/styles" -Method Post -Body $body -ContentType "application/json" | Out-Null

$networkId = (Invoke-RestMethod -Uri "$baseUrl/networks")[0]
# Apply the style
Invoke-RestMethod -Uri "$baseUrl/apply/styles/ReferenceImageStyle/$networkId" -Method Get | Out-Null

# Apply a standard, unclustered force-directed layout
Invoke-RestMethod -Uri "$baseUrl/commands/layout/force-directed" -Method Post -Body '{}' -ContentType "application/json" | Out-Null
Start-Sleep -Seconds 3

# Export image
Invoke-RestMethod -Uri "$baseUrl/networks/$networkId/views/first.png?h=2500" -OutFile "d:\Zika_wetlab\meta-analysis\plots\11_Cytoscape_Reference_Style.png"
