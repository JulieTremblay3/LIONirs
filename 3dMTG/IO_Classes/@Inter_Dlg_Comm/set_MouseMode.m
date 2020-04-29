function oInterDlgComm = set_MouseMode( oInterDlgComm, strMode )

    ItemNo = getfield( oInterDlgComm.MouseMode.Subscripts, strMode );
    oInterDlgComm.MouseData.Mode = ItemNo;
