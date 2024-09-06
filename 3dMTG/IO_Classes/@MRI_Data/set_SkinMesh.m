function oMRI = set_SkinMesh( oMRI, VertexBuffer, IndexBuffer  )
    oMRI.Surfaces.Skin.VertexBuffer = VertexBuffer;
    oMRI.Surfaces.Skin.IndexBuffer  = IndexBuffer;