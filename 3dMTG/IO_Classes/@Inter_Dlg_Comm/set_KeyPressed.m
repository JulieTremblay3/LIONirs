function oObj = set_KeyPressed( oObj, KeyLowerCase )
    oObj.KeyboardData.LastKey = KeyLowerCase;
    oObj.KeyboardData.LKClock = clock;