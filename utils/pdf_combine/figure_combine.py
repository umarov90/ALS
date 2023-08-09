import joblib
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
from params import Params
import anndata as ad
import PyPDF2

p = Params()

# Set the figure size (adjust as needed)
fig_width = 8
fig_height = 6

# Initialize the PDF
pdf_pages = []

sub_info = "day"

for info_val in ["D4", "D7", "D14"]:
    original_path = f"{p.folder}figures/scvelo_day_{info_val}.png"
    # Load the images using PIL
    original_image = Image.open(original_path)

    # Create a new figure
    fig, axs = plt.subplots(1, 1, figsize=(fig_width, fig_height))

    # Plot the images side by side
    axs.imshow(original_image)
    axs.axis('off')

    # Add a common title for the row (equal to ID)
    fig.suptitle(f'Day: {info_val}', fontsize=12)

    # Save the figure as a PDF page
    buf = BytesIO()
    plt.savefig(buf, format='pdf')
    buf.seek(0)
    pdf_pages.append(buf)

    # Close the current figure to avoid memory leaks
    plt.close()

# Combine all PDF pages into a single PDF file
with open('../output.pdf', 'wb') as output_file:
    pdf_writer = PyPDF2.PdfWriter()
    for page in pdf_pages:
        pdf_reader = PyPDF2.PdfReader(page)
        pdf_writer.add_page(pdf_reader.pages[0])
    pdf_writer.write(output_file)
