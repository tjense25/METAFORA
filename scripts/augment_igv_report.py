import sys
import os
import re
import base64
from pathlib import Path
from pdf2image import convert_from_path

############################################
# READ IN COMMAND LINE SETTINGS
############################################
PDF_DIR = sys.argv[1]
PRIOT_OUTS = sys.argv[2]
INPUT_HTML = sys.argv[3]
OUTPUT_HTML = sys.argv[4]
############################################

def pdf_to_jpg_base64(pdf_path):
    """
    Convert first page of PDF to JPG and return base64 string.
    """
    pages = convert_from_path(pdf_path, dpi=120)
    img = pages[0]

    jpg_path = pdf_path.replace(".pdf",".jpg")
    img.save(jpg_path, "JPEG", quality=60)

    with open(jpg_path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode("utf-8")

    os.remove(jpg_path)
    return f"data:image/jpeg;base64,{encoded}"


def extract_filename_from_region(region):
    """
    Expect regions like:
        chr11:61025584-61028394
    Returns:
        chr11_61025584_61028394.methylation_outlier_plots.pdf
    """
    match = re.search(r"(chr[\w]+):(\d+)-(\d+)", region)
    if match:
        chrom, start, end = match.groups()
        return f"{chrom}_{start}_{end}.methylation_outlier_plots.pdf"
    return None


def build_plot_dictionary(prioritized_outliers, outlier_dir):
    plot_dict={}
    for i,coord in enumerate(prioritized_outliers):
        pdf = os.path.join(outlier_dir, extract_filename_from_region(coord))
        if os.path.exists(pdf):
            print(f"Processing {pdf}")
            img_string = pdf_to_jpg_base64(pdf)
            plot_dict[i] = img_string
    return plot_dict

def generate_injected_block(plot_dict):
    """
    Generate HTML+JS block to inject at bottom of HTML.
    """

    # Build JS dictionary string
    js_dict_entries = []
    for index, img in plot_dict.items():
        js_dict_entries.append(f'"{index}": "{img}"')

    js_dict_string = ",\n".join(js_dict_entries)

    injected_block = f"""
<!-- ===================== -->
<!-- OUTLIER PLOT PANEL -->
<!-- ===================== -->
<div id="outlierPlotContainer" style="
    padding: 12px 16px;
    border-top: 1px solid #ccc;
    border-bottom: 1px solid #ccc;
    background: white;
">
    <div style="font-size: 12px; font-weight: 600; margin: 0 0 8px 0;">Outlier plot</div>
    <img id="outlierPlot"
         style="width:95%; max-width:1000; max-height:400; height:auto; display:block;"
         src="">
</div>
<script>
const plotMap = {{
{js_dict_string}
}};
function updateOutlierPlot(uniqueId) {{
    if (plotMap[uniqueId]) {{
        document.getElementById("outlierPlot").src = plotMap[uniqueId];
    }}
}}
updateOutlierPlot("0");
</script>
"""
    return injected_block


def inject_into_html(input_html, output_html, injected_block):
    with open(input_html, "r", encoding="utf-8") as f:
        html = f.read()

    if "</body>" not in html:
        raise RuntimeError("No </body> tag found in HTML.")

    new_html = html.replace("</body>", injected_block + "\n</body>")
    new_html = new_html.replace("igvBrowser.loadSession({", "console.log(uniqueId);\nupdateOutlierPlot(uniqueId)\nigvBrowser.loadSession({")
    layout_script = r"""
<script>
window.addEventListener("load", function() {
  const container = document.getElementById("container");
  const table = document.getElementById("tableContainer");
  const igv = document.getElementById("igvContainer");
  const plot = document.getElementById("outlierPlotContainer");

  if (!container || !table || !igv || !plot) {
    console.warn("Could not reorder layout:", {container, table, igv, plot});
    return;
  }

  // Ensure vertical layout
  container.style.display = "flex";
  container.style.flexDirection = "column";

  // Remove plot from wherever it currently is (likely end of body)
  if (plot.parentNode) plot.parentNode.removeChild(plot);

  // Insert plot right after the table, before IGV
  // (table is inside container; igv is also inside container)
  if (table.nextSibling) {
    container.insertBefore(plot, table.nextSibling);
  } else {
    container.appendChild(plot);
  }

  // Make sure IGV stays after the plot
  if (igv.parentNode !== container) {
    if (igv.parentNode) igv.parentNode.removeChild(igv);
    container.appendChild(igv);
  }

  // Optional: prevent plot from shrinking weirdly
  plot.style.flex = "0 0 auto";
  table.style.flex = "0 0 auto";
  igv.style.flex = "1 1 auto";
});
</script>
"""
    new_html = new_html.replace("</body>", layout_script + "\n</body>")

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(new_html)

def main():
    prioritized_outliers=[row.strip().split()[-1] for row in open(PRIOT_OUTS,'r')][1:]
    plot_dict = build_plot_dictionary(prioritized_outliers,PDF_DIR)
    injected_block = generate_injected_block(plot_dict)
    inject_into_html(INPUT_HTML, OUTPUT_HTML, injected_block)
    print(f"\nDone. Output written to: {OUTPUT_HTML}")

if __name__ == "__main__":
    main()
