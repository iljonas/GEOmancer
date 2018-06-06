<!DOCTYPE html>
<html>
<body>

<style type='text/css'>
	.sampleheader{border:1px solid black;color:#000;font-weight:1000;text-align:center;height:50px;}
    .samplerow{border:1px solid black;color:#000;font-weight:500;;text-align:center;height:30px;width:107px}
    .mainheader{border:1px solid black;color:#000;font-weight:1000;text-align:left;height:30px;}
    .mainrow{border:1px solid black;color:#000;font-weight:500;;text-align:left;height:30px;width:235px;padding-left:10px}
</style>

<table style="border-spacing:1;">
	<tr>
		<label><td class="mainheader">Study ID: </td>
		<td class="mainrow"><input type="text" id="study id" size="8"></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">Platform ID: </td>
		<td class="mainrow"><input type="text" id="platform id" size="8"></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">Outside study controls: </td>
		<td class="mainrow"><input type="text" id="study controls" size="30"></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">Difference Calculation: </td>
		<td class="mainrow"><select id="diff calc">
			<option value="Subtraction" SELECTED>Subtraction
			<option value="Log2">Log Base 2
			</select></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">Comparison Calculation: </td>
		<td class="mainrow"><select id="comp calc">
			<option value="Mean" SELECTED>Mean
			<option value="Sum">Sum
			</select></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">In-column Delimiter: </td>
		<td class="mainrow"><input type="text" id="delim" size="2"></td>
		</label>
	</tr>
	<tr>
		<label><td class="mainheader">Enter number of Samples: </td>
		<td class="mainrow"><input type="number" id="sample count" size="4"></td>
		</label>
	</tr>
</table>
<p><span id='display'></span></p>
<button onclick="ValidateTable();">Show Samples</button>

<script>

function generateTable(rowNum){
	document.write('<table style="border-spacing:0;" id="Table2">')
	document.write('	<thead>')
	document.write('		<tr>')
	document.write('      <td class="sampleheader">Sample ID</td>')
	document.write('			<td class="sampleheader">Group</td>')
	document.write('			<td class="sampleheader">Species</td>')
	document.write('			<td class="sampleheader">Strain</td>')
	document.write('			<td class="sampleheader">Treatment</td>')
	document.write('			<td class="sampleheader">Dosage</td>')
	document.write('			<td class="sampleheader">Exposure Time</td>')
	document.write('			<td class="sampleheader">Growth Phase</td>')
	document.write('			<td class="sampleheader">Other</td>')
	document.write('		</tr>')
	document.write('	</thead>')
	document.write('	<tbody>')
	for (i=0;i<rowNum;i++){
		document.write('<tr>')
        document.write('	<td class="samplerow"><input type="text" id="sample id" size="10" pattern="GSM([0-9]{4,5})"></td>')
        document.write('	<td class="samplerow"><select id="group" style="width: 100px">')
		document.write('		<option value="Control">Control')
		document.write('		<option value="Test">Test')
		document.write('		<option value="Neither" SELECTED>Neither')
		document.write('	</select></td>')
		document.write('	<td class="samplerow"><input type="text" id="species" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="strain" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="treatment" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="dosage" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="exposure time" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="growth phase" size="10"></td>')
		document.write('	<td class="samplerow"><input type="text" id="other" size="10"></td>')
		document.write('</tr>')  
	}
    document.write('	</tbody>')
	document.write('</table>')
	document.write('<p id="sample errors">HELLO</p>')
	document.write('<button onclick="exportTableToCSV(\'output.csv\');">Export</button>')
}

function downloadCSV(csv, filename) {
    var csvFile;
    var downloadLink;

    // CSV file
    csvFile = new Blob([csv], {type: "text/csv"});

    // Download link
    downloadLink = document.createElement("a");

    // File name
    downloadLink.download = filename;

    // Create a link to the file
    downloadLink.href = window.URL.createObjectURL(csvFile);

    // Hide download link
    downloadLink.style.display = "none";

    // Add the link to DOM
    document.body.appendChild(downloadLink);

    // Click downljavascript:;oad link
    downloadLink.click();
}

function exportTableToCSV(filename) {
		var csv = [];
    var rows = document.querySelectorAll("tbody tr");
    var tbl = [];
    var cellVal;
    var regex = new RegExp("GSM([0-9]{4-7})");
    for (var i = 0; i < rows.length; i++) {
		var row = [], cols = rows[i].querySelectorAll("td, select, input");
        for (var j = 0; j < cols.length; j++) {
       		if(j == 1 && (!cols[j].value.startsWith("GSM") || cols[j].value == "")) {
				document.getElementById("sample errors").innerHTML = "Sample " + (i + 1) + " needs to be formatted as GSM12345";
        		return;
	       	} else if(j == 1 && cols[j].value.startsWith("GSM")) {
            	document.getElementById("sample errors").innerHTML = "";
            }
			if(cols[j].value == "") {
				cellVal = "NA";
			} else if(cols[j].value == undefined && cols[j].innerText.length > 0) {
				cellVal = cols[j].innerText;
			} else if(cols[j].value != "") {
				cellVal = cols[j].value;
			}
			if(cellVal != undefined) {
				row.push(cellVal);
        	}
        }
		csv.push(row.join());        
    }

    // Download CSV file
   //downloadCSV(csv.join("\n"), filename);
}

function ValidateTable(){
	var errorText = "";
	if(document.getElementById("study id").value == "" || document.getElementById("platform id").value == "" || document.getElementById("sample count").value == "") {
		errorText = "Please enter the Study ID, Platform ID, and Sample Count";
	} else if(!document.getElementById("study id").value.match("GSE([0-9]{4,5})")) {
		errorText = "Please ensure Study ID is in this format: GSE12345";
	} else if(!document.getElementById("platform id").value.match("GPL([0-9]{4,5})")) {
		errorText = "Please ensure Platform ID is in this format: GPL12345";
	} else if(document.getElementById("delim").value.length > 1) {
		errorText = "Please enter a one character delimiter or leave blank";
	}
	if (errorText != ""){
		document.getElementById("display").innerHTML = errorText;
	} else {
		generateTable(document.getElementById("sample count").value);
	}
}             
</script>


</body>
</html>