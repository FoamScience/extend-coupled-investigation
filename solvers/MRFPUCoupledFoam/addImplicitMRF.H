
    if (addMRF)
    {
        // Add MRF contribution to block diagonal
        const tensorField& TUIn = tTU().internalField();
        CoeffField<vector4>::squareTypeField& DD = UpEqn.diag().asSquare();
        const scalarField& V = mesh.V().field();

        forAll (mrfZones, mrfZoneI)
        {
            const labelList& curZoneCells = mrfZones[mrfZoneI].zone();

            forAll (curZoneCells, zcI)
            {
                label cellI = curZoneCells[zcI];
                const scalar& cellV = V[cellI];
                const tensor& cellTU = TUIn[cellI];
                CoeffField<vector4>::squareType& cellDD = DD[cellI];

                // Add MRF tensor to velocity block (3x3 upper-left)
                cellDD(0, 0) += cellV*cellTU.xx();
                cellDD(0, 1) += cellV*cellTU.xy();
                cellDD(0, 2) += cellV*cellTU.xz();

                cellDD(1, 0) += cellV*cellTU.yx();
                cellDD(1, 1) += cellV*cellTU.yy();
                cellDD(1, 2) += cellV*cellTU.yz();

                cellDD(2, 0) += cellV*cellTU.zx();
                cellDD(2, 1) += cellV*cellTU.zy();
                cellDD(2, 2) += cellV*cellTU.zz();
            }
        }
    }
