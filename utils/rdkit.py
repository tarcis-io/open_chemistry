def create_compound_images(compounds : list, image_size : int, directory : str):
    """
    Creates images of chemical compounds using RDKit library.

    Parameters:
        - compounds  (list) : List of compound dictionaries with keys 'CID' and 'CanonicalSMILES'.
        - image_size (int)  : Size of images.
        - directory  (str)  : Directory where the images will be saved.

    Returns:
        - None

    Dependencies:
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
