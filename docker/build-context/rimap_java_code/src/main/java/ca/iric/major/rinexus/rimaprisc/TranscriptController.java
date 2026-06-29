/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class TranscriptController {

    @GetMapping(value = "/visualize", produces = "text/html")
    public String getVisualization() {
        // Example data
        int transcriptLength = 4107;
        int[] interactionPositions = {3719, 1619, 3108};

        // HTML structure with Canvas
        StringBuilder html = new StringBuilder();
        html.append("<!DOCTYPE html>");
        html.append("<html lang='en'>");
        html.append("<head>");
        html.append("<meta charset='UTF-8'>");
        html.append("<meta name='viewport' content='width=device-width, initial-scale=1.0'>");
        html.append("<title>Transcript Visualization</title>");
        html.append("</head>");
        html.append("<body>");
        html.append("<h1>miRNA Interaction Sites on Transcript</h1>");
        html.append("<canvas id='transcriptCanvas' width='800' height='200' style='border:1px solid black;'></canvas>");
        html.append("<script>");
        html.append("const canvas = document.getElementById('transcriptCanvas');");
        html.append("const ctx = canvas.getContext('2d');");
        html.append("const transcriptLength = ").append(transcriptLength).append(";"); // Dynamic length
        html.append("const interactionPositions = ").append(java.util.Arrays.toString(interactionPositions)).append(";"); // Dynamic positions
        html.append("const margin = 50;");
        html.append("const width = canvas.width - 2 * margin;");
        html.append("const startY = canvas.height / 2;");

        // Draw transcript line
        html.append("ctx.beginPath();");
        html.append("ctx.moveTo(margin, startY);");
        html.append("ctx.lineTo(canvas.width - margin, startY);");
        html.append("ctx.strokeStyle = 'black';");
        html.append("ctx.lineWidth = 2;");
        html.append("ctx.stroke();");

        // Draw interaction markers
        html.append("interactionPositions.forEach(position => {");
        html.append("    const x = margin + (position / transcriptLength) * width;");
        html.append("    ctx.beginPath();");
        html.append("    ctx.moveTo(x, startY - 10);");
        html.append("    ctx.lineTo(x, startY + 10);");
        html.append("    ctx.strokeStyle = 'red';");
        html.append("    ctx.lineWidth = 2;");
        html.append("    ctx.stroke();");
        html.append("    ctx.fillStyle = 'black';");
        html.append("    ctx.font = '12px Arial';");
        html.append("    ctx.fillText(position, x - 10, startY + 20);");
        html.append("});");
        html.append("</script>");
        html.append("</body>");
        html.append("</html>");

        return html.toString();
    }
}
