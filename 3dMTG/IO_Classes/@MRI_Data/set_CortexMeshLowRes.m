function oMRI = set_CortexMeshLowRes( oMRI, VertexBuffer, IndexBuffer )
    oMRI.Surfaces.CortexLowRes.VertexBuffer = VertexBuffer;
    oMRI.Surfaces.CortexLowRes.IndexBuffer  = IndexBuffer;
    