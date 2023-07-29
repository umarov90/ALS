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

original_path = f"{p.folder}figures/scvelo__original_ALL.png"
meta_path = f"{p.folder}figures/scvelo__meta_ALL.png"

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
fig.suptitle("scvelo", fontsize=12)

plt.savefig(p.folder + "scvelo_pair.pdf")

plt.close()
