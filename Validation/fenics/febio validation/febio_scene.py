import numpy as np 
import meshio, tempfile, subprocess




def febio_scene(parameters, elements): 
    # Create the FEBio spec file

    c = []
    c.append(f'<?xml version="1.0" encoding="UTF-8"?>')
    c.append(f'<febio_spec version="3.0">')
    c.append(f'  <Module type="solid"/>')
    c.append(f'  <Control>')
    c.append(f'    <time_steps>{parameters["nsteps"]}</time_steps>')
    c.append(f'    <step_size>1</step_size>')
    c.append(f'    <solver>')
    c.append(f'    	<max_refs>{parameters["nnewtonsteps"]}</max_refs>')
    c.append(f'    	<max_ups>{parameters["nBFGSsteps"]}</max_ups>')
    c.append(f'    	<diverge_reform>1</diverge_reform>')
    c.append(f'    	<reform_each_time_step>1</reform_each_time_step>')
    c.append(f'    	<dtol>{parameters["displacement_tolerance"]}</dtol>')
    c.append(f'    	<etol>0.01</etol>')
    c.append(f'    	<rtol>{parameters["residual_tolerance"]}</rtol>')
    c.append(f'    	<lstol>0</lstol>')
    c.append(f'    	<min_residual>1e-15</min_residual>')
    c.append(f'    	<qnmethod>BFGS</qnmethod>')
    c.append(f'    	<rhoi>-2</rhoi>')
    c.append(f'    </solver>')
    c.append(f'  </Control>')
    c.append(f'  <Material>')
    c.append(f'    <material id="1" name="{parameters["febio_material"]}" type="{parameters["febio_material"]}">')
    if(parameters['febio_material'] != "Mooney-Rivlin"): 
        c.append(f'      <E>{parameters["young_modulus"]}</E>')
        c.append(f'      <v>{parameters["poisson_ratio"]}</v>')
    else: 
        c.append(f'     <c1>{parameters["c01"]}</c1>')
        c.append(f'     <c2>{parameters["c10"]}</c2>')
        c.append(f'     <k>{parameters["k"]}</k>')
    
    c.append(f'    </material>')
    c.append(f'  </Material>')

    c.append(f'  <Geometry>')
    c.append(f'    <Nodes name="nodes">')
    for id, n in enumerate(parameters["mesh"].points):
        c.append(f'      <node id="{id+1}">{",".join([str(i) for i in n])}</node>')
    c.append(f'    </Nodes>')
    c.append(f'    <NodeSet name="base_nodes">')
    base_nodes = np.unique(np.concatenate([elements[parameters["bas_id"]][type] for type in elements[parameters["bas_id"]].keys()])).astype(int)
    for n in base_nodes:
        c.append(f'      <node id="{n+1}"/>')
    c.append(f'    </NodeSet>')

    c.append(f'    <NodeSet name="top_nodes">')
    top_nodes = np.unique(np.concatenate([elements[parameters["top_id"]][type] for type in elements[parameters["top_id"]].keys()])).astype(int)
    for n in top_nodes:
        c.append(f'      <node id="{n+1}"/>')
    c.append(f'    </NodeSet>')
   
    for type, indices in elements[parameters["vol_id"]].items():
        c.append(f'    <Elements type="{type}" mat="1">')
        for id, t in enumerate(indices[0].tolist()):
            c.append(f'      <elem id="{id+1}">{",".join([str(i+1) for i in t])}</elem>')
        c.append(f'    </Elements>')
    c.append(f'    <Surface name="base">')
    id = 0
    for type, indices in elements[parameters["bas_id"]].items():
        for t in indices.tolist():
            c.append(f'      <{type} id="{id+1}">{",".join([str(i+1) for i in t])}</{type}>')
            id += 1
    c.append(f'    </Surface>')
    c.append(f'    <Surface name="top">')
    id = 0
    for type, indices in elements[parameters["top_id"]].items():
        for t in indices.tolist():
            c.append(f'      <{type} id="{id+1}">{",".join([str(i+1) for i in t])}</{type}>')
            id += 1
    c.append(f'    </Surface>')
    c.append(f'  </Geometry>')

    c.append(f'  <Boundary>')
    c.append(f'  <bc name="FixedDisplacement01" type="fix" node_set="base_nodes">')
    c.append(f'  		<dofs>x,y,z</dofs>')
    c.append(f'  </bc>')
    c.append(f'  </Boundary>')

    c.append(f'  <Loads>')
    c.append(f'      <nodal_load name="ForceLoad01" type="nodal_load" node_set="top_nodes">')
    c.append(f'      		<dof>y</dof>')
    c.append(f'      		<scale lc="1">{parameters["traction"][1]}</scale>')
    c.append(f'      	</nodal_load>')
    c.append(f'  </Loads>')

    c.append(f'  <LoadData>')
    c.append(f'  	<load_controller id="1" type="loadcurve">')
    c.append(f'  		<interpolate>SMOOTH</interpolate>')
    c.append(f'  		<points>')
    c.append(f'  			<point>0,0</point>')
    c.append(f'  			<point>1,1</point>')
    c.append(f'  		</points>')
    c.append(f'  	</load_controller>')
    c.append(f'  </LoadData>')

    c.append(f'  <Output>')
    c.append(f'    <logfile file="{parameters["displacement_filename"]}">')
    c.append(f'      <node_data data="x;y;z"></node_data>')
    c.append(f'    </logfile>')
    c.append(f'  </Output>')

    c.append(f'</febio_spec>')


    with tempfile.NamedTemporaryFile(suffix=".feb") as f:
        feb_filename = f.name
    with open(feb_filename, 'w') as f:
        f.write('\n'.join(c) + '\n')

    print('Running febio...', end='')
    args = [
        "-i",
        feb_filename
    ]
    p = subprocess.Popen(
        ['febio3'] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    while True:
        line = p.stdout.readline()
        if not line:
            break
        # line = line.decode("utf-8")
        # print(line, end = '')

    p.communicate()
    assert p.returncode == 0, f"febio exited with error (return code {p.returncode})."
    print('Done')
    with open(parameters["displacement_filename"], 'r') as f:
        lines = f.readlines()
        last_record_position = 0
        for id, l in zip(range(len(lines)),lines):
            if l.startswith("Data Record"):
                last_record_position = id

        start = last_record_position+5
        end = start + len(parameters["mesh"].points)
        displacement_lines = lines[start:end]

        displacements = np.zeros((len(displacement_lines), 3))
        for id, l in zip(range(len(displacement_lines)),displacement_lines):
            l = l[:-1].split(' ')[1:]
            displacements[id] = np.array([l[0], l[1], l[2]]).astype(float)

    assert len(displacements) == len(parameters["mesh"].points)
    return displacements

# Sort element indices per tag and type
def sort_elements(meshio_to_febio_element_types, parameters):
    elements = {}
    for block, tags in zip(parameters["mesh"].cells, parameters["mesh"].cell_data['gmsh:physical']):
        type, indices = block
        type = meshio_to_febio_element_types[type]
        for t in np.unique(tags):
            if t not in elements:
                elements[t] = {}
            mask = np.array(np.ma.masked_not_equal(tags, t).mask)
            if len(indices[~mask]):
                if type not in elements[t]:
                    elements[t][type] = indices[~mask]
                else:
                    elements[t][type] = np.concatenate((indices[~mask], elements[t][type]), axis=0)

    for tag, types in elements.items():
        if len(types) > 1:
            print(f'Tag "{tag}" has "')
            for type, indices in types.items():
                print(f'\t{len(indices)} elements of type "{type}"')
        else:
            type, indices = list(types.items())[0]
            print(f'Tag "{tag}" has {len(indices)} elements of type "{type}"')

    print("Mesh has {} nodes".format(len(parameters["mesh"].points)))
    return elements

