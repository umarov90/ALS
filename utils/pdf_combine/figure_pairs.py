import joblib
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
from params import Params
import PyPDF2

p = Params()

# Set the figure size (adjust as needed)
fig_width = 10
fig_height = 5

# Initialize the PDF
pdf_pages = []

ids = joblib.load(f"{p.folder}ids.p")

for ID in ids:
    original_path = f"{p.folder}figures/scvelo_original_{ID}.png"
    meta_path = f"{p.folder}figures/scvelo_meta_{ID}.png"

    # Load the images using PIL
    original_image = Image.open(original_path)
    meta_image = Image.open(meta_path)

    # Create a new figure
    fig, axs = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    # Plot the images side by side
    axs[0].imshow(original_image)
    axs[0].axis('off')
    axs[0].set_title('Original', fontsize=10)

    axs[1].imshow(meta_image)
    axs[1].axis('off')
    axs[1].set_title('Meta', fontsize=10)

    # Add a common title for the row (equal to ID)
    fig.suptitle(f'ID: {ID}', fontsize=12)

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
