from ezdxf import readfile

def dxf_polygon(file_path):
    """
    Extracts vertices from DXF entities of LINE, LWPOLYLINE, or POLYLINE types in the specified DXF file.

    Parameters:
        file_path (str): The path to the DXF file.

    Returns:
        list: A list of tuples representing the extracted vertices.

    Raises:
        ValueError: If the lines are not connected.
        ValueError: If there are multiple POLYLINE entities and the function does not support them.
        ValueError: If the entities in the file are neither all lines nor all polylines.
    """
    doc = readfile(file_path)
    ents = list(doc.modelspace().query('LINE LWPOLYLINE POLYLINE'))
    ents_type = [e.DXFTYPE for e in ents]
    if all(e == 'LINE' for e in ents_type):
        vertices = []
        ents.append(ents[0])
        for i in range(len(ents)-1):
            p0 = (ents[i].dxf.end.x, ents[i].dxf.end.y)
            p1 = (ents[i+1].dxf.start.x, ents[i+1].dxf.start.y)
            if p0 == p1:
                vertices.append(p0)
            else:
                raise ValueError(f"Lines are not connected! Check LINE#{ents[i].dxf.handle} and LINE#{ents[i+1].dxf.handle}.")
            
    elif all(e == 'LWPOLYLINE' for e in ents_type):
        if len(ents) > 1:
            raise ValueError("Does not support multiple POLYLINE.")

        vertices = [(q[0],q[1]) for q in ents[0].get_points()]
       
    elif all(e == 'POLYLINE' for e in ents_type):
        if len(ents) > 1:
            raise ValueError("Does not support multiple POLYLINE.")
        
        vertices = [(q.dxf.location.x,q.dxf.location.y) for q in ents[0].vertices]
        
    else:
        raise ValueError("The entities in the file are not all lines or all polylines.")   
    return vertices