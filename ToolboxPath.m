pwd
path = pwd;
Mac=1;
if Mac==1; sep='/'; else; sep='\'; end
path_sim=strcat(path,sep,'Simulations');
path_app=strcat(path,sep,'HCP_Application');
save('ToolboxPath','path','Mac')
save(strcat(path_sim,sep,'ToolboxPath'),'path','Mac')
save(strcat(path_app,sep,'ToolboxPath'),'path','Mac')