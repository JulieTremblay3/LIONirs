function Ret = get_CurKeyPressed( oObj )

    e = etime(clock, oObj.KeyboardData.LKClock);
    
    if( e < 0.6 )
        Ret = oObj.KeyboardData.LastKey;
    else
        Ret = 0;
    end
