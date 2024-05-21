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
    image_file = image.filename.rsplit('.', 1)[0]

    for transposition in transpositions:

        image_with_transposition = image.transpose(transposition)
        image_with_transposition.save('{}_with_transposition_{}.png'.format(image_file, transposition))
