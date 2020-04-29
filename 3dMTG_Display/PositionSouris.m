%Fonction permettant de connaitre la position du pointeur de la souris
function [pointerI, in] = PositionSouris(figureH, itemH)    
    in=0;
    resolution = get(0, 'ScreenSize');
    pointerS = get(0, 'PointerLocation')./resolution(3:4);
    posFig = get(figureH, 'Position');
    pointerF = (pointerS - posFig(1:2))./posFig(3:4);
    posItem = get(itemH, 'Position');
    pointerA = (pointerF - posItem(1:2))./posItem(3:4);
    if (sum(pointerA < 0) == 0) & (sum(pointerA > 1) == 0)
        in= 1;
    end
    pointerI = pointerA;