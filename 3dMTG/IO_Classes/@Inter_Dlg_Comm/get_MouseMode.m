function bActive = get_MouseModeActive( oInterDlgComm, strMode )

    ItemNo = getfield( oInterDlgComm.MouseMode.Subscripts, strMode );
    return (ItemNo==oInterDlgComm.MouseData.Mode);
