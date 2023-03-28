def Heterostructure(Materials,sz=5,h=3.5,thetas=[0.0*sc.DEGREE]):
    '''
    **CURRENTLY ONLY WORKS WITH TWO LAYERS

    INPUT:  Materials:   List of files to be in the heterostructure
                         ***add file type checker for now cif or vasp

            h:           Separation between layers
                         ** will add support to accept lists

            sz:          maximum number of sub lattices along 'a' and 'b'

            thetas:      angles between layers. list object or scalar

    OUTPUT: Final_ase:   Heterostructure ase
            strain:      strain value
            Final_sc:    Heterostructure supercell-core

    - Need to test different shifts along (a,b), the supercell may account for several possible shifts depending on its size

    *** Need to adjust substrate to be at the origin
    '''
    C_sup = 0.0
    for j, material in enumerate(Materials):
        if material.endswith('.vasp'):
            mat_ase = ase.io.vasp.read_vasp(material)
        elif material.endswith('.cif'):
            mat_ase = ase.io.vasp.read_cif(material)
        else:
            mat_ase = ase.io.vasp.read_vasp(material)

        C_sup += mat_ase.get_cell_lengths_and_angles()[2]
    layers_ase = []
    layers_sc  = []
    n = len(Materials)
    for j, material in enumerate(Materials):

        mat_ase = ase.io.vasp.read_vasp(material)
        atom_pos = mat_ase.get_scaled_positions()
        U = mat_ase.get_cell()

        mat_sc = sc.lattice()
        #Create Matrix for Unit cell based on a,b,c,alpha,beta,gamma

        # Default units are Angstrom coordinates
        mat_sc.set_vectors(U[0],U[1],U[2])

        d = 0.0
        s = 0.0

        Ci =  mat_ase.cell.cellpar()[2]
        if j == 0:
            h_sub = np.max(mat_ase.positions[:,2])+ 0.001
        if  j>0:
            s  = (h_sub)/Ci

        # Change it so substrate layer starts at zero
        for pos in atom_pos:
            if pos[2] >=0.5 or pos[2] <0.:
                if pos[2] >= 0.5:
                    diff = 1 - pos[2]
                elif pos[2] < 0.:
                    diff = 0 - pos[2]
                if diff > d:
                    d = diff

        ele = mat_ase.get_chemical_symbols()
        for i, pos in enumerate(atom_pos):
            if pos[2] >=0.5:
                pos[2] = (pos[2] - 1 + d + s)*Ci/C_sup
                atom_pos[i,2] = pos[2]
                mat_sc.add_atom(ele[i],[pos[0],pos[1],pos[2]], unit=sc.Unit.Crystal)
            elif pos[2] <0.5:
                pos[2] = (pos[2] + d + s)*Ci/C_sup
                atom_pos[i,2] = pos[2]
                mat_sc.add_atom(ele[i],[pos[0],pos[1],pos[2]], unit=sc.Unit.Crystal)
        mat_ase.set_scaled_positions(atom_pos)

        layers_ase.append(mat_ase)
        layers_sc.append(mat_sc)
    supercell = sc.heterostructure()
    supercell.set_substrate(layers_sc[0])
    for i in range(1,len(layers_ase)):
        supercell.add_layer(layers_sc[i])

    # n-1 layers because the substrate isn't counted for supercell core

    tmp = supercell.opt(max_el=sz, thetas=thetas)
    final = supercell.calc(M=tmp.M(), thetas=tmp.thetas())

    #Gives the Heterostructure Unit Cell
    N = len(final.superlattice().atoms())


    ## Create Final Heterostructure in ASE
    ele_sup = []
    pos_sup = []
    for atom_i in final.superlattice().atoms():
        ele_sup.append(atom_i.element)
        pos_sup.append(atom_i.pos)

    U = final.superlattice().vectors()
    Final_ase = Atoms(cell = [(U[0][0],U[0][1],U[0][2]),(U[1][0],U[1][1],U[1][2]),(U[2][0],U[2][1],U[2][2])],
                     positions = pos_sup,
                     symbols   = ele_sup )
    index = 0
    smol  = 999

    for i, pos in enumerate(Final_ase.positions):
        if pos[2] > h_sub and index == 0:
            index = i
            smol = np.min(Final_ase.positions[index:,2])
        if pos[2] > h_sub:
            pos[2] = pos[2] - smol + h + h_sub
    # WE Need the height of each layer plus the number of layers to properly adjust more than two layers

    strain = final.max_strain()
    Final_sc = final.superlattice()
    return([Final_ase,strain,Final_sc])
