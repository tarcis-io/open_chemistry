def create_images_with_transpositions(image_file : str):
    """
    Creates images with transpositions from the provided image file.

    Parameters:
        - image_file (str) : Image file in which the transposed images will be created.

    Returns:
        - None

    Dependencies:
        - pillow
    """

    from PIL import Image

    transpositions = [
        Image.FLIP_LEFT_RIGHT,
        Image.FLIP_TOP_BOTTOM,
        Image.ROTATE_90,
        Image.ROTATE_180,
        Image.ROTATE_270,
        Image.TRANSPOSE,
        Image.TRANSVERSE
    ]

    image      = Image.open(image_file)
    image_file = image_file.rsplit('.', 1)[0]

    for transposition in transpositions:

        image_with_transposition = image.transpose(transposition)
        image_with_transposition.save('{}_with_transposition_{}.png'.format(image_file, transposition))


def plot_images_with_transpositions(image_file : str, title : str):
    """
    Plots the provided image file and its transposed images.

    Parameters:
        - image_file (str) : Image file that has transposed images.
        - title      (str) : Title of the image.

    Returns:
        - None

    Dependencies:
        - math
        - matplotlib
        - pillow
    """

    import math
    import matplotlib.pyplot
    from PIL import Image

    transpositions = [
        Image.FLIP_LEFT_RIGHT,
        Image.FLIP_TOP_BOTTOM,
        Image.ROTATE_90,
        Image.ROTATE_180,
        Image.ROTATE_270,
        Image.TRANSPOSE,
        Image.TRANSVERSE
    ]

    files = [
        image_file
    ]

    image_file = image_file.rsplit('.', 1)[0]

    for transposition in transpositions:

        files.append('{}_with_transposition_{}.png'.format(image_file, transposition))

    fig_cols = 4
    fig_rows = math.ceil(len(files) / fig_cols)

    fig, axs = matplotlib.pyplot.subplots(fig_rows, fig_cols, figsize = (16, 8), constrained_layout = True)
    fig.suptitle(title, fontsize = 16)

    axs = axs.flatten()
    for ax in axs: ax.axis('off')

    for index, file in enumerate(files):

        ax = axs[index]
        ax.imshow(matplotlib.pyplot.imread(file))
