function [VertexBuffer, IndexBuffer] = get_CortexMeshHiRes( oMRI )
    VertexBuffer = oMRI.Surfaces.CortexHiRes.VertexBuffer;
    IndexBuffer = oMRI.Surfaces.CortexHiRes.IndexBuffer;