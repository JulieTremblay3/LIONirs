function [VertexBuffer, IndexBuffer] = get_SkinMesh( oMRI )
    VertexBuffer = oMRI.Surfaces.Skin.VertexBuffer;
    IndexBuffer = oMRI.Surfaces.Skin.IndexBuffer;