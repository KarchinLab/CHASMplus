---
#
# Use the widgets beneath and the content will be
# inserted automagically in the webpage. To make
# this work, you have to use › layout: frontpage
#
layout: frontpage
#header:
#  image_fullwidth: chasmplus.png
#
# Use the call for action to show a button on the frontpage
#
# To make internal links, just use a permalink like this
# url: /getting-started/
#
# To style the button in different colors, use no value
# to use the main color or success, alert or secondary.
# To change colors see sass/_01_settings_colors.scss
#
permalink: /index.html
---
<div style="text-align:center"><h1>Explore driver mutations</h1></div>

<div id="tabulator-controls" class="table-controls">
  <div class="medium-8" style="display:block;margin-left:auto;margin-right:auto;">
    <input name="name" type="text" placeholder="Search by gene name" class="form-control">
  </div>
  <div class="medium-2 columns">
    <div style="float:left;">Cancer type:</div>
    <select name="ctype" placeholder="Select a cancer type">
      <option value="All">All types</option>
      <option value="THCA">THCA</option>
      <option value="BRCA">BRCA</option>
      <option value="LGG">LGG</option>
      <option value="UCEC">UCEC</option>
      <option value="GBM">GBM</option>
      <option value="LIHC">LIHC</option>
      <option value="STAD">STAD</option>
      <option value="PRAD">PRAD</option>
      <option value="BLCA">BLCA</option>
      <option value="OV">OV</option>
      <option value="LUAD">LUAD</option>
      <option value="UVM">UVM</option>
      <option value="PAAD">PAAD</option>
      <option value="LUSC">LUSC</option>
      <option value="COAD">COAD</option>
      <option value="UCS">UCS</option>
      <option value="TGCT">TGCT</option>
      <option value="READ">READ</option>
      <option value="HNSC">HNSC</option>
      <option value="KIRP">KIRP</option>
      <option value="KIRC">KIRC</option>
      <option value="MESO">MESO</option>
      <option value="ESCA">ESCA</option>
      <option value="CESC">CESC</option>
      <option value="LAML">LAML</option>
      <option value="CHOL">CHOL</option>
      <option value="SARC">SARC</option>
      <option value="DLBC">DLBC</option>
      <option value="THYM">THYM</option>
      <option value="KICH">KICH</option>
      <option value="ACC">ACC</option>
      <option value="PCPG">PCPG</option>
    </select>
  </div>
  <div class="medium-6 columns" style="display:block;margin-left:auto;margin-right:auto;">
    Columns:<br>
    <input type="checkbox" name="Hugo_Symbol" checked> Gene
    <input type="checkbox" name="Transcript_ID" checked> Transcript
    <input type="checkbox" name="HGVSp_Short" checked> Mutation
    <input type="checkbox" name="gwCHASMplus score" checked> Score
    <input type="checkbox" name="cancer type with highest prevalence" checked> Cancer type
    <br>
    <input type="checkbox" name="number of mutations (highest cancer type)" checked> Number of mutations
    <input type="checkbox" name="frequency category (highest cancer type)" checked> Frequency
    <input type="checkbox" name="url" checked> Detailed information
  </div>
  <div class="medium-4 columns">
    Download:<br>
    <button type="button" name="csv-download">CSV</button>
    <button type="button" name="xlsx-download">Excel</button>
  </div>
</div>
<div id="example-table" class="medium-12 columns"></div>

<!--
<br>
<br>
<p><strong>Interactive viewers:</strong></p>
<ul class="side-nav">
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_205745">ACC: Adrenocortical carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_205936">BLCA: Bladder Urothelial Carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_205945">BRCA: Breast invasive carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210009">CESC: Cervical squamous cell carcinoma and endocervical adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210000">CHOL: Cholangiocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210104">COAD: Colon adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210049">DLBC: Lymphoid Neoplasm Diffuse Large B-cell Lymphoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210127">ESCA: Esophageal carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210116">GBM: Glioblastoma multiforme</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210138">HNSC: Head and Neck squamous cell carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210213">KICH: Kidney Chromophobe</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210149">KIRC: Kidney renal clear cell carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210247">KIRP: Kidney renal papillary cell carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210234">LAML: Acute Myeloid Leukemia</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210225">LGG: Brain Lower Grade Glioma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210300">LIHC: Liver hepatocellular carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210313">LUAD: Lung adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210339">LUSC: Lung squamous cell carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210326">MESO: Mesothelioma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210352">OV: Ovarian serous cystadenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210401">PAAD: Pancreatic adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210424">PANCAN: Pan-cancer (multiple cancer types)</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210416">PCPG: Pheochromocytoma and Paraganglioma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210447">PRAD: Prostate adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210439">READ: Rectum adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210432">SARC: Sarcoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210453">SKCM: Skin Cutaneous Melanoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210527">STAD: Stomach adenocarcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210517">TGCT: Testicular Germ Cell Tumors</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210507">THCA: Thyroid carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210500">THYM: Thymoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210552">UCEC: Uterine Corpus Endometrial Carcinoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210542">UCS: Uterine Carcinosarcoma</a></li>
<li><a href="http://www.cravat.us/CRAVAT/job_detail.html?job_id=collintokheim_20180815_210535">UVM: Uveal Melanoma</a></li>
</ul>
-->
