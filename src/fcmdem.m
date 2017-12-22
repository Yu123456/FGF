function fcmdem(action)

global FcmFigh FcmFigTitle FcmAxish FcmCenter FcmU OldDataID;
if nargin==0
    action='initialize';
end
if strcmp(action,'initialize')
    FcmFigTitle='2-D Fuzzy C-Means Clustering';
    FcmFigH=findobj(0,'Name',FcmFigTitle);
    if isempty(FcmFigH)
        eval([mfilename,'(''set_gui'')'])
        set(findobj(FcmFigH,'Units','pixels'),'Units','normal')
        set(findobj(FcmFigH,'Interrupt','off'),'Interrupt','on')
    else
        refresh(FcmFigH)
    end
elseif strcmp(action,'set_gui')
    FcmFigH=figure('Name',FcmFigTitle,'NumberTitle','off','DockControls','off','Resize','off')
end


end

