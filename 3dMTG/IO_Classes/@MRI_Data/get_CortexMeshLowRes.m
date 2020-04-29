function [VertexBuffer, IndexBuffer] = get_CortexMeshLowRes( oMRI )
    VertexBuffer = oMRI.Surfaces.CortexLowRes.VertexBuffer;
    IndexBuffer = oMRI.Surfaces.CortexLowRes.IndexBuffer;