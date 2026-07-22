$baseUrl = "http://localhost:1234/v1"

# Create Style
$style = @{
    title = "ZikaPublicationStyle2"
    defaults = @(
        @{ visualProperty = "NODE_SHAPE"; value = "ELLIPSE" },
        @{ visualProperty = "NODE_SIZE"; value = 45 },
        @{ visualProperty = "NODE_LABEL_COLOR"; value = "#000000" },
        @{ visualProperty = "NODE_LABEL_FONT_SIZE"; value = 14 },
        @{ visualProperty = "NODE_LABEL_POSITION"; value = "C,C,c,0,0" },
        @{ visualProperty = "NETWORK_BACKGROUND_PAINT"; value = "#FFFFFF" },
        @{ visualProperty = "EDGE_WIDTH"; value = 2 },
        @{ visualProperty = "EDGE_TRANSPARENCY"; value = 150 }
    )
    mappings = @(
        @{
            mappingType = "passthrough"
            mappingColumn = "label"
            mappingColumnType = "String"
            visualProperty = "NODE_LABEL"
        },
        @{
            mappingType = "passthrough"
            mappingColumn = "color"
            mappingColumnType = "String"
            visualProperty = "NODE_FILL_COLOR"
        }
    )
}

$body = $style | ConvertTo-Json -Depth 10
Invoke-RestMethod -Uri "$baseUrl/styles" -Method Post -Body $body -ContentType "application/json" | Out-Null

# Apply Style
$networkId = (Invoke-RestMethod -Uri "$baseUrl/networks")[0]
Invoke-RestMethod -Uri "$baseUrl/apply/styles/ZikaPublicationStyle2/$networkId" -Method Get | Out-Null

# Export Image
$exportBody = @{
    options = "PNG"
    OutputFile = "d:\Zika_wetlab\meta-analysis\plots\Cytoscape_Publication_Network.png"
} | ConvertTo-Json

Invoke-RestMethod -Uri "$baseUrl/commands/view/export" -Method Post -Body $exportBody -ContentType "application/json" | Out-Null
