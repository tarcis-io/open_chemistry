def create_compound_images(compounds : list, directory : str, image_size : int):
    """
    Creates images of chemical compounds using the RDKit library.

    Parameters:
        - compounds  (list) : List of the compound dictionaries with the keys 'CID' and 'CanonicalSMILES'.
        - directory  (str)  : Directory where the images will be saved.
        - image_size (int)  : Size of the images.

    Returns:
        - None

    Dependencies:
        - os
        - rdkit
    """

    import os
    from rdkit import Chem

    for compound in compounds:

        CID             = compound['CID']
        CanonicalSMILES = compound['CanonicalSMILES']

        image_file = f'{ CID }_RDKit.png'
        image_file = os.path.join(directory, image_file)

        molecule = Chem.MolFromSmiles(CanonicalSMILES)
        Chem.Draw.MolToFile(molecule, image_file, size = (image_size, image_size))


def plot_compound_images(compounds : list, directory : str, title : str):
    """
    Plots the images of the chemical compounds created by the RDKit library.

    Parameters:
        - compounds (list) : List of the compound dictionaries with the keys 'CID' and 'Title'.
        - directory (str)  : Directory where the images are saved.
        - title     (str)  : Title of the image.

    Returns:
        - None

    Dependencies:
        - math
        - matplotlib
        - os
    """

    import math
    import matplotlib.pyplot
    import os

    fig_cols = 5
    fig_rows = math.ceil(len(compounds) / fig_cols)

    fig, axs = matplotlib.pyplot.subplots(fig_rows, fig_cols, figsize = (16, 16), constrained_layout = True)
    fig.suptitle(title, fontsize = 16)

    axs = axs.flatten()
    for ax in axs: ax.axis('off')

    for index, compound in enumerate(compounds):

        CID   = compound['CID']
        Title = compound['Title']

        image_file = f'{ CID }_RDKit.png'

        for root, dirs, files in os.walk(directory):

            if image_file in files:

                image_file = os.path.join(root, image_file)
                break

        ax = axs[index]
        ax.set_title(Title)
        ax.imshow(matplotlib.pyplot.imread(image_file))
