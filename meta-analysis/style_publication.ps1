$baseUrl = "http://localhost:1234/v1"

# Create a true Publication-Ready Style
$style = @{
    title = "NaturePublicationStyle"
    defaults = @(
        @{ visualProperty = "NODE_SHAPE"; value = "ELLIPSE" },
        @{ visualProperty = "NODE_BORDER_WIDTH"; value = 0.0 },
        @{ visualProperty = "NODE_LABEL_COLOR"; value = "#222222" },
        @{ visualProperty = "NODE_LABEL_FONT_SIZE"; value = 14 },
        @{ visualProperty = "NODE_LABEL_POSITION"; value = "S,C,c,0.00,5.00" },
        @{ visualProperty = "NETWORK_BACKGROUND_PAINT"; value = "#FFFFFF" },
        @{ visualProperty = "EDGE_WIDTH"; value = 1.0 },
        @{ visualProperty = "EDGE_STROKE_UNSELECTED_PAINT"; value = "#CCCCCC" },
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
            mappingType = "continuous"
            mappingColumn = "degree"
            mappingColumnType = "Double"
            visualProperty = "NODE_SIZE"
            points = @(
                @{ value = 1.0; lesser = "20.0"; equal = "20.0"; greater = "20.0" },
                @{ value = 90.0; lesser = "100.0"; equal = "100.0"; greater = "100.0" }
            )
        },
        @{
            mappingType = "continuous"
            mappingColumn = "degree"
            mappingColumnType = "Double"
            visualProperty = "NODE_FILL_COLOR"
            points = @(
                @{ value = 1.0; lesser = "#FEE0D2"; equal = "#FEE0D2"; greater = "#FEE0D2" },
                @{ value = 40.0; lesser = "#FC9272"; equal = "#FC9272"; greater = "#FC9272" },
                @{ value = 90.0; lesser = "#DE2D26"; equal = "#DE2D26"; greater = "#DE2D26" }
            )
        },
        @{
            mappingType = "continuous"
            mappingColumn = "degree"
            mappingColumnType = "Double"
            visualProperty = "NODE_LABEL_FONT_SIZE"
            points = @(
                @{ value = 1.0; lesser = "1"; equal = "1"; greater = "1" },
                @{ value = 40.0; lesser = "12"; equal = "12"; greater = "12" },
                @{ value = 90.0; lesser = "24"; equal = "24"; greater = "24" }
            )
        }
    )
}

$body = $style | ConvertTo-Json -Depth 10
Invoke-RestMethod -Uri "$baseUrl/styles" -Method Post -Body $body -ContentType "application/json" | Out-Null

$networkId = (Invoke-RestMethod -Uri "$baseUrl/networks")[0]
# Apply the style
Invoke-RestMethod -Uri "$baseUrl/apply/styles/NaturePublicationStyle/$networkId" -Method Get | Out-Null

# Apply a publication-ready layout (Force-Directed usually looks very nice with proper weighting, or Kamada-Kawai)
Invoke-RestMethod -Uri "$baseUrl/commands/layout/kamada-kawai" -Method Post -Body '{}' -ContentType "application/json" | Out-Null
Start-Sleep -Seconds 3

# Export image
Invoke-RestMethod -Uri "$baseUrl/networks/$networkId/views/first.png?h=2500" -OutFile "d:\Zika_wetlab\meta-analysis\plots\10_Cytoscape_Nature_Style.png"
