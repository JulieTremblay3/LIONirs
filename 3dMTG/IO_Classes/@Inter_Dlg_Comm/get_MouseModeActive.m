function bActive = get_MouseModeActive( oInterDlgComm, strMode )

    ItemNo = getfield( oInterDlgComm.MouseMode.Subscripts, strMode );
    bActive = (ItemNo==oInterDlgComm.MouseData.Mode);
